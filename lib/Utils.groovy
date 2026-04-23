//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

import org.yaml.snakeyaml.Yaml

class Utils {

    //
    // When running with -profile conda, warn if channels have not been set-up appropriately
    //
    public static void checkCondaChannels(log) {
        Yaml parser = new Yaml()
        def channels = []
        try {
            def config = parser.load("conda config --show channels".execute().text)
            channels = config.channels
        } catch(NullPointerException | IOException e) {
            log.warn "Could not verify conda channel configuration."
            return
        }

        // Check that all channels are present
        // This channel list is ordered by required channel priority.
        def required_channels_in_order = ['conda-forge', 'bioconda', 'defaults']
        def channels_missing = ((required_channels_in_order as Set) - (channels as Set)) as Boolean

        // Check that they are in the right order
        def channel_priority_violation = false
        def n = required_channels_in_order.size()
        for (int i = 0; i < n - 1; i++) {
            channel_priority_violation |= !(channels.indexOf(required_channels_in_order[i]) < channels.indexOf(required_channels_in_order[i+1]))
        }

        if (channels_missing | channel_priority_violation) {
            log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  There is a problem with your Conda configuration!\n\n" +
                "  You will need to set-up the conda-forge and bioconda channels correctly.\n" +
                "  Please refer to https://bioconda.github.io/\n" +
                "  The observed channel order is \n" +
                "  ${channels}\n" +
                "  but the following channel order is required:\n" +
                "  ${required_channels_in_order}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        }
    }

    public static remove_lanes_from_meta(meta_map, lane_keys = ['num_lanes', 'lane', 'lanes', 'id', 'read_group', 'size', 'data_type']) {
        def new_meta = meta_map - meta_map.subMap(lane_keys) + [id: meta_map.sample]
        return new_meta
    }

    // public static robustly_test_if_empty(value) {
    //     def truthy_val = (
    //         (!value) ||
    //         (value == null) || 
    //         (value instanceof Collection && value.empty) || 
    //         (value instanceof String && value.replaceAll(/["']/, "").trim() == "") ||
    //         (value instanceof File && value.isEmpty())
    //     )
    //     return truthy_val
    // }
    static boolean robustly_test_if_empty(value) {
        if (value == null) return true

        // Collections / Maps / arrays: emptiness is definitive — do not fall through
        // to the File check, since stringifying "[Path(/x)]" does not name a real file.
        if (value instanceof Collection) return value.isEmpty()
        if (value instanceof Map) return value.isEmpty()
        if (value instanceof Object[]) return value.length == 0

        // CharSequence covers String and GString; treat whitespace/quote-only as empty
        if (value instanceof CharSequence) {
            return value.toString().replaceAll(/["']/, "").trim() == ""
        }

        // File / Path: empty iff missing, zero-length file, or empty directory
        if (value instanceof File) {
            if (!value.exists()) return true
            if (value.isFile() && value.length() == 0) return true
            if (value.isDirectory() && value.list().length == 0) return true
            return false
        }
        if (value instanceof java.nio.file.Path) {
            def f = value.toFile()
            if (!f.exists()) return true
            if (f.isFile() && f.length() == 0) return true
            if (f.isDirectory() && f.list().length == 0) return true
            return false
        }

        // Groovy-truthy fallback: 0, false, etc. → empty
        if (!value) return true

        // Iterables without a specific type above — try to realize to a list
        try {
            def val = value.toList()
            if (val instanceof List) return val.isEmpty()
        } catch (Exception ignored) { }

        // Last resort: treat as a path-like string
        def f = new File(value.toString())
        if (!f.exists()) return true
        if (f.isFile() && f.length() == 0) return true
        if (f.isDirectory() && f.list().length == 0) return true

        return false
    }
}
