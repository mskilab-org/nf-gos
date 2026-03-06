// memcap.c
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dlfcn.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdarg.h>

static long get_mem_limit_kb() {
    char *mem_limit_str = getenv("MEMLIMIT_MB");
    if (!mem_limit_str) return 0;
    return atol(mem_limit_str) * 1024;
}

static int create_fake_meminfo() {
    long mem_limit_kb = get_mem_limit_kb();
    if (mem_limit_kb == 0) return -1;
    
    char template[] = "/tmp/meminfo_XXXXXX";
    int fd = mkstemp(template);
    if (fd < 0) return -1;
    
    char buf[2048];
    int len = snprintf(buf, sizeof(buf),
        "MemTotal:       %ld kB\n"
        "MemFree:        %ld kB\n"
        "MemAvailable:   %ld kB\n"
        "Buffers:        0 kB\n"
        "Cached:         0 kB\n"
        "SwapCached:     0 kB\n"
        "Active:         0 kB\n"
        "Inactive:       0 kB\n"
        "SwapTotal:      0 kB\n"
        "SwapFree:       0 kB\n",
        mem_limit_kb,
        mem_limit_kb * 9 / 10,
        mem_limit_kb * 9 / 10);
    
    write(fd, buf, len);
    unlink(template);
    lseek(fd, 0, SEEK_SET);
    return fd;
}

int open(const char *pathname, int flags, ...) {
    static int (*real_open)(const char*, int, ...) = NULL;
    if (!real_open) {
        real_open = dlsym(RTLD_NEXT, "open");
    }
    
    if (pathname && strcmp(pathname, "/proc/meminfo") == 0) {
        int fd = create_fake_meminfo();
        if (fd >= 0) return fd;
    }
    
    mode_t mode = 0;
    if (flags & O_CREAT) {
        va_list args;
        va_start(args, flags);
        mode = va_arg(args, int);
        va_end(args);
        return real_open(pathname, flags, mode);
    }
    return real_open(pathname, flags);
}

int open64(const char *pathname, int flags, ...) {
    static int (*real_open64)(const char*, int, ...) = NULL;
    if (!real_open64) {
        real_open64 = dlsym(RTLD_NEXT, "open64");
    }
    
    if (pathname && strcmp(pathname, "/proc/meminfo") == 0) {
        int fd = create_fake_meminfo();
        if (fd >= 0) return fd;
    }
    
    mode_t mode = 0;
    if (flags & O_CREAT) {
        va_list args;
        va_start(args, flags);
        mode = va_arg(args, int);
        va_end(args);
        return real_open64(pathname, flags, mode);
    }
    return real_open64(pathname, flags);
}

int openat(int dirfd, const char *pathname, int flags, ...) {
    static int (*real_openat)(int, const char*, int, ...) = NULL;
    if (!real_openat) {
        real_openat = dlsym(RTLD_NEXT, "openat");
    }
    
    if (pathname && strcmp(pathname, "/proc/meminfo") == 0) {
        int fd = create_fake_meminfo();
        if (fd >= 0) return fd;
    }
    
    mode_t mode = 0;
    if (flags & O_CREAT) {
        va_list args;
        va_start(args, flags);
        mode = va_arg(args, int);
        va_end(args);
        return real_openat(dirfd, pathname, flags, mode);
    }
    return real_openat(dirfd, pathname, flags);
}

FILE *fopen(const char *pathname, const char *mode) {
    static FILE *(*real_fopen)(const char*, const char*) = NULL;
    if (!real_fopen) {
        real_fopen = dlsym(RTLD_NEXT, "fopen");
    }
    
    if (pathname && strcmp(pathname, "/proc/meminfo") == 0) {
        int fd = create_fake_meminfo();
        if (fd >= 0) return fdopen(fd, "r");
    }
    
    return real_fopen(pathname, mode);
}

long sysconf(int name) {
    static long (*real_sysconf)(int) = NULL;
    if (!real_sysconf) {
        real_sysconf = dlsym(RTLD_NEXT, "sysconf");
    }
    
    long mem_limit_kb = get_mem_limit_kb();
    if (mem_limit_kb == 0) {
        return real_sysconf(name);
    }
    
    long page_size = real_sysconf(_SC_PAGESIZE);
    
    switch(name) {
        case _SC_PHYS_PAGES:
        case _SC_AVPHYS_PAGES:
            return (mem_limit_kb * 1024) / page_size;
        default:
            return real_sysconf(name);
    }
}