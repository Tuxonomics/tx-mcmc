// utils.h


#include <stdio.h>
#include <execinfo.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
// #include <math.h>
#include <cmath>


#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795028841971
#endif


#ifndef M_SQRT1_2
	#define M_SQRT1_2 0.707106781186547524401
#endif


#if !defined(TUXONOMICS_BASIC_TYPES)
    #define TUXONOMICS_BASIC_TYPES

    typedef uint8_t  u8;
    typedef uint16_t u16;
    typedef uint32_t u32;
    typedef uint64_t u64;

    typedef int8_t  i8;
    typedef int16_t i16;
    typedef int32_t i32;
    typedef int64_t i64;

    typedef float  f32;
    typedef double f64;

    typedef i8  b8;
    typedef i32 b32;
#endif

    typedef float  s;
    typedef double d;

#if !defined(TUXONOMICS_UTILITIES)

    #define TUXONOMICS_UTILITIES

    #if defined(_MSC_VER)
        #if _MSC_VER < 1300
            #define DEBUG_TRAP() __asm int 3
        #else
            #define DEBUG_TRAP() __debugbreak()
        #endif
    #else
        #define DEBUG_TRAP() __builtin_trap()
    #endif

    #if !defined(TEST)
    #if !defined(RELEASE) && !defined(ASSERTS)
        #define ASSERT_MSG_VA(cond, msg, ...) do { \
            if (!(cond)) { \
            assertHandler(__FILE__, (i32)__LINE__, msg, __VA_ARGS__); \
            DEBUG_TRAP(); \
            } \
            } while(0)

        #define ASSERT_MSG(cond, msg) ASSERT_MSG_VA(cond, msg, 0)

        #define ASSERT(cond) ASSERT_MSG_VA(cond, 0, 0)
        #define PANIC(msg) ASSERT_MSG_VA(0, msg, 0)
        #define UNIMPLEMENTED() ASSERT_MSG_VA(0, "unimplemented", 0);
    #else
        #define ASSERT_MSG_VA(cond, msg, ...)
        #define ASSERT_MSG(cond, msg)
        #define ASSERT(cond)
        #define PANIC(msg)
        #define UNIMPLEMENTED()
    #endif
    #endif


    #if !defined(INLINE)
        #if defined(_MSC_VER)
            #if _MSC_VER < 1300
                #define INLINE
            #else
                #define INLINE __forceinline
            #endif
        #else
            #define INLINE __attribute__ ((__always_inline__))
        #endif
    #endif


    #if !defined(_Threadlocal)
        #if defined(_MSC_VER)
            #define _Threadlocal __declspec( thread )
        #else
            #define _Threadlocal __thread
        #endif
    #endif


    void Backtrace() {
    #define BACKTRACE_MAX_STACK_DEPTH 50
    #if SYSTEM_POSIX
        void* callstack[BACKTRACE_MAX_STACK_DEPTH];
        int i, frames = backtrace(callstack, BACKTRACE_MAX_STACK_DEPTH);
        char** strs = backtrace_symbols(callstack, frames);
        for (i = 0; i < frames; ++i) {
            fprintf(stderr, "%s\n", strs[i]);
        }
        free(strs);
    #elif SYSTEM_WINDOWS
        UNIMPLEMENTED();
    #endif
    }

    void assertHandler(char const *file, i32 line, char const *msg, ...) {
        va_list args;
        va_start(args, msg);
        Backtrace();

        if (msg) {
            fprintf(stderr, "Assert failure: %s:%d: ", file, line);
            vfprintf(stderr, msg, args);
            fprintf(stderr, "\n");
        } else {
            fprintf(stderr, "Assert failure: %s:%d\n", file, line);
        }
        va_end(args);
    }

#endif


#ifndef MIN
    #define MIN(x, y) ((x) <= (y) ? (x) : (y))
    #define MAX(x, y) ((x) >= (y) ? (x) : (y))
#endif


#ifndef TUXONOMICS_ALLOCATOR
    #define TUXONOMICS_ALLOCATOR

    typedef enum AllocType {
        AT_Alloc,
        AT_Calloc,
        AT_Realloc,
        AT_Free,
        AT_FreeAll,
    } AllocType;


    #define ALLOC_FUNC(name) void *name(void *payload, enum AllocType alType, size_t count, size_t size, void *old)
    typedef void *allocFunc(void *payload, enum AllocType alType, size_t count, size_t size, void *old);


    typedef struct Allocator Allocator;
    struct Allocator {
        allocFunc *func;
        void *payload;
    };


    void *Alloc(Allocator al, size_t count) {
        return al.func(al.payload, AT_Alloc, count, 0, NULL);
    }

    void *Calloc(Allocator al, size_t count, size_t size) {
        return al.func(al.payload, AT_Calloc, count, size, NULL);
    }

    void *Free(Allocator al, void* ptr) {
        if (ptr)
            al.func(al.payload, AT_Free, 0, 0, ptr);
        return NULL;
    }

    void *FreeAll(Allocator al) {
        al.func(al.payload, AT_FreeAll, 0, 0, NULL);
        return NULL;
    }

    void *Realloc(Allocator al, void *ptr, size_t size, size_t oldsize) {
        return al.func(al.payload, AT_Realloc, size, oldsize, ptr);
    }


    void *checkedCalloc(size_t num_elems, size_t elem_size) {
        void *ptr = calloc(num_elems, elem_size);
        if (!ptr) {
            perror("calloc failed");
            exit(1);
        }
        return ptr;
    }

    void *checkedRealloc(void *ptr, size_t num_bytes) {
        ptr = realloc(ptr, num_bytes);
        if (!ptr) {
            perror("realloc failed");
            exit(1);
        }
        return ptr;
    }

    void *checkedMalloc(size_t num_bytes) {
        void *ptr = malloc(num_bytes);
        if (!ptr) {
            perror("malloc failed");
            exit(1);
        }
        return ptr;
    }

    void *heapAllocFunc(void *payload, enum AllocType alType, size_t count, size_t size, void *old) {
        switch (alType) {
            case AT_Alloc:
                return checkedMalloc(count);
            case AT_Calloc:
                return checkedCalloc(count, size);
            case AT_Free:
            case AT_FreeAll: {
                free(old);
                return NULL;
            }
            case AT_Realloc:
                return checkedRealloc(old, count);
        }
        return NULL;
    }

    Allocator _ALLOCATOR_DEFAULT = { .func = heapAllocFunc, .payload = 0 };


    typedef struct Arena {
        Allocator allocator;
        u8  *raw;
        u64 cap;
        u64 len;
    } Arena;

    void *arenaAllocFunc(void *payload, enum AllocType alType, size_t count, size_t size, void *old) {
        Arena *arena = (Arena *) payload;

        switch (alType) {
            case AT_Alloc: {
                if (arena->len + count > arena->cap) {
                    return NULL;
                }
                u8 *ptr = &arena->raw[arena->len];
                arena->len += count;
                return ptr;
            }
            case AT_Calloc: {
                u8 * ptr = (u8 *) arenaAllocFunc( payload, AT_Alloc, count * size, 0, old );
                memset( ptr, 0, (count * size) );
                return ptr;
            }
            case AT_Free:
            case AT_FreeAll: {
                arena->len = 0;
                break;
            }
            case AT_Realloc: {
                break;
            }
        }

        return NULL;
    }

    Allocator ArenaAllocatorMake(Arena *arena) {
        Allocator al;
        al.func    = arenaAllocFunc;
        al.payload = arena;
        return al;
    }

    void ArenaInit(Arena *arena, Allocator al, u64 size) {
        arena->allocator = al;
        arena->raw = (u8 *) Alloc(al, size);
        arena->cap = size;
        arena->len = 0;
    }

    void ArenaDefaultInit(Arena *arena, u64 size) {
        ArenaInit(arena, _ALLOCATOR_DEFAULT, size);
    }

    void ArenaDestroy(Arena *arena) {
        if ( arena->raw )
            Free( arena->allocator, arena->raw );
    }

    void ArenaDestroyResize( Arena *arena, u64 size ) {
        ArenaDestroy( arena );

        arena->raw = (u8 *) Alloc( arena->allocator, size );
        arena->cap = size;
        arena->len = 0;
    }

    void ArenaInitAndAllocator( Allocator alOnArena, Arena *arena, Allocator *al, u64 size )
    {
        ArenaInit( arena, alOnArena, size );
        *al = ArenaAllocatorMake( arena );
    }

    void ArenaAllocatorCheck( Arena *arena, Allocator *al, u64 size )
    {
        if ( arena->cap < size ) {
            if ( arena->raw ) {
                ArenaDestroyResize( arena, size );
            }
            else {
                ArenaInitAndAllocator( _ALLOCATOR_DEFAULT, arena, al, size );
            }
        }
    }
#endif

#if TEST
void test_arena()
{
    u32 aSize = 2;

    Arena arena;
    ArenaDefaultInit( &arena, aSize * sizeof( u32 ) );

    Allocator arenaAllocator = ArenaAllocatorMake( &arena );

    u32 *a = (u32 *) Alloc(  arenaAllocator, sizeof( u32 ) );
    u32 *b = (u32 *) Calloc( arenaAllocator, 1, sizeof( u32 ) );

    ASSERT( ((u64) b - (u64) a) == 4 );

    u32 bSize = 3;

    ArenaDestroyResize( &arena, bSize * sizeof( u32 ) );

    a = (u32 *) Alloc(  arenaAllocator, sizeof( u32 ) );
    b = (u32 *) Calloc( arenaAllocator, 1, sizeof( u32 ) );

    u32 *c = (u32 *) Alloc(  arenaAllocator, sizeof( u32 ) );

    ASSERT( ((u64) c - (u64) b) == 4 );

    ArenaDestroy( &arena );
}
#endif


INLINE
b32 Equal( f64 a, f64 b, f64 eps )
{
	if ( fabs(a - b) <= eps ) {
		return 1;
	}
	return 0;
}


// from http://www.gingerbill.org/article/2015/08/19/defer-in-cpp/
template <typename F>
struct privDefer {
    F f;
    privDefer(F f) : f(f) {}
    ~privDefer() { f(); }
};

template <typename F>
privDefer<F> defer_func(F f) {
    return privDefer<F>(f);
}

#define DEFER_1(x, y) x##y
#define DEFER_2(x, y) DEFER_1(x, y)
#define DEFER_3(x)    DEFER_2(x, __COUNTER__)
#define defer(code)   auto DEFER_3(_defer_) = defer_func([&](){code;})


