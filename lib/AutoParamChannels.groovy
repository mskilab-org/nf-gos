// lib/AutoParamChannels.groovy

@Singleton
class AutoParamChannels {
    private Map channels = [:]
    private Map channelTypes = [:]  // Track whether file or value channel
    private boolean initialized = false
    
    /**
     * Auto-initialize all params as appropriate channel types
     */
    def get(String name) {
        if (!initialized) {
            initializeFromParams()
        }
        
        if (!channels.containsKey(name)) {
            error "Parameter '${name}' not found. Available: ${channels.keySet()}"
        }
        
        return channels[name]
    }
    
    /**
     * Initialize all channels from params automatically
     */
    private void initializeFromParams() {
        params.each { key, value ->
            if (value != null && value != '') {
                channels[key] = createAppropriateChannel(key, value)
            }
        }
        initialized = true
    }
    
    /**
     * Determine if value is a file path/URI and create appropriate channel
     */
    private createAppropriateChannel(String key, value) {
        // Skip special Nextflow params
        if (isSystemParam(key)) {
            channelTypes[key] = 'system'
            return Channel.value(value)
        }
        
        // Handle lists/collections
        if (value instanceof List || value instanceof Collection) {
            return createChannelFromList(key, value)
        }
        
        // Check if it's a file path or URI
        if (isFilePath(value)) {
            channelTypes[key] = 'file'
            return createFileChannel(value)
        }
        
        // Default to value channel
        channelTypes[key] = 'value'
        return Channel.value(value)
    }
    
    /**
     * Check if a value looks like a file path or URI
     */
    private boolean isFilePath(value) {
        if (!(value instanceof String)) {
            return false
        }
        
        String str = value.toString().trim()
        
        // Check for URI pattern (anything://...)
        if (str =~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\//) {
            return true
        }
        
        // Also check for specific known cloud storage URIs (for clarity/documentation)
        if (str.startsWith('s3://') || 
            str.startsWith('gs://') || 
            str.startsWith('az://') ||
            str.startsWith('gcs://') ||
            str.startsWith('http://') ||
            str.startsWith('https://') ||
            str.startsWith('ftp://')) {
            return true
        }
        
        // Check for absolute paths
        if (str.startsWith('/') || str.startsWith('~/')) {
            return true
        }
        
        // Check for relative paths with file extensions
        if (str.contains('/') || hasFileExtension(str)) {
            return true
        }
        
        // Check if file exists locally (catches relative paths)
        try {
            def file = new File(str)
            if (file.exists() || file.parent != null) {
                return true
            }
        } catch (Exception e) {
            // Not a valid file path
        }
        
        return false
    }
    
    /**
     * Check if string has a file extension
     */
    private boolean hasFileExtension(String str) {
        // First check if there's ANY extension (anything after a dot in the last part)
        def lastPart = str.split('/')[-1]  // Get filename part
        if (lastPart.contains('.') && !lastPart.startsWith('.')) {
            // Has at least one extension
            return true
        }
        
        return false
    }
    
    /**
     * Create file channel with proper handling for cloud storage
     */
    private createFileChannel(filePath) {
        try {
            // Nextflow's fromPath handles cloud URIs automatically
            return Channel.fromPath(filePath, checkIfExists: true).collect()
        } catch (Exception e) {
            // If file doesn't exist but looks like a path, create empty channel
            println "Warning: File not found for param: ${filePath}"
            return Channel.empty()
        }
    }
    
    /**
     * Handle list parameters
     */
    private createChannelFromList(String key, list) {
        // Check if all items look like files
        if (list.every { isFilePath(it) }) {
            channelTypes[key] = 'file_list'
            return Channel.fromPath(list).collect()
        }
        
        // Otherwise treat as value list
        channelTypes[key] = 'value_list'
        return Channel.value(list)
    }
    
    /**
     * Check if parameter is a Nextflow system param
     */
    private boolean isSystemParam(String key) {
        def systemParams = [
            'help', 'resume', 'profile', 'work-dir', 'config',
            'revision', 'latest', 'with-tower', 'with-report',
            'with-timeline', 'with-trace', 'with-dag',
            'max_memory', 'max_cpus', 'max_time'
        ]
        
        return systemParams.any { key.startsWith(it) || key.contains('-') }
    }
    
    /**
     * Get type of channel for a parameter
     */
    String getType(String name) {
        if (!initialized) {
            initializeFromParams()
        }
        return channelTypes[name] ?: 'unknown'
    }
    
    /**
     * Get all file channels
     */
    Map getFileChannels() {
        if (!initialized) {
            initializeFromParams()
        }
        return channels.findAll { k, v -> channelTypes[k] == 'file' }
    }
    
    /**
     * Get all value channels
     */
    Map getValueChannels() {
        if (!initialized) {
            initializeFromParams()
        }
        return channels.findAll { k, v -> channelTypes[k] == 'value' }
    }
    
    /**
     * Check if parameter exists and is not empty
     */
    boolean has(String name) {
        if (!initialized) {
            initializeFromParams()
        }
        return channels.containsKey(name) && channels[name] != Channel.empty()
    }
    
    /**
     * Direct property access
     */
    def propertyMissing(String name) {
        return get(name)
    }
}

// Convenience accessor
class PC {
    static def get(String name) {
        return AutoParamChannels.instance.get(name)
    }
    
    static boolean has(String name) {
        return AutoParamChannels.instance.has(name)
    }
    
    static String type(String name) {
        return AutoParamChannels.instance.getType(name)
    }
}