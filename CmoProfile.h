#ifndef CMO_PROFILE_H
#define CMO_PROFILE_H

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <x86intrin.h>
// #include <threads.h>
#include <pthread.h>

typedef uint32_t bool32_t;

typedef struct ProfileEntry
{
    const char* func;
    uint64_t count;
    bool32_t exit;
} ProfileEntry;

typedef struct ProfileScratch
{
    int threadId;
    int end;
    int maxCount;
    ProfileEntry* entries;
} ProfileScratch;

// extern _Thread_local ProfileScratch threadProf;
extern _Thread_local int CmoProfileId;
extern ProfileScratch* profiler;

void cmo_allocate_prof_array(int numThreads, int entries);
void cmo_free_prof_array(int numThreads);
void cmo_init_prof(ProfileScratch* scr, int threadId, int entries);
void cmo_free_prof(ProfileScratch* scr);
void cmo_prof_push_entry(ProfileScratch* scr, const char* name, 
                         uint64_t count, bool32_t exit);
void cmo_prof_start_timed_region(ProfileScratch* scr, const char* name);
void cmo_prof_end_timed_region(ProfileScratch* scr, const char* name);
void cmo_prof_print_file(ProfileScratch* scr, const char* fname);

#ifndef CMO_NO_PROF
#define CMO_PROF_FUNC_START() cmo_prof_start_timed_region(&profiler[CmoProfileId], __func__)
#define CMO_PROF_FUNC_END() cmo_prof_end_timed_region(&profiler[CmoProfileId], __func__)
#define CMO_PROF_REGION_START(reg) cmo_prof_start_timed_region(&profiler[CmoProfileId], reg)
#define CMO_PROF_REGION_END(reg) cmo_prof_end_timed_region(&profiler[CmoProfileId], reg)
#else
#define CMO_PROF_FUNC_START()
#define CMO_PROF_FUNC_END()
#define CMO_PROF_REGION_START(reg)
#define CMO_PROF_REGION_END(reg)
#endif


#ifdef CMO_PROFILE_IMPL
#include <stdio.h>

// https://github.com/cmuratori/meow_hash/issues/31
#ifndef NOT_PIC_SAFE
static uint64_t cmo_measure_start(void)
{
    uint32_t a, b;
    __asm__ __volatile__(
        "push %%rbx\n\t"
        "cpuid\n\t"
        "rdtsc\n\t"
        "mov %%eax, %0\n\t"
        "mov %%edx, %1\n\t"
        "pop %%rbx" : "=r"(a), "=r"(b) : : "rax", "rcx", "rdx");
    return ((uint64_t)b << 32) | a;
}

static uint64_t cmo_measure_end(void)
{
    uint32_t a, b;
    __asm__ __volatile__(
        "push %%rbx\n\t"
        "rdtscp\n\t"
        "mov %%eax, %0\n\t"
        "mov %%edx, %1\n\t"
        "cpuid\n\t"
        "pop %%rbx" : "=r"(a), "=r"(b) : : "rax", "rcx", "rdx");
    return ((uint64_t)b << 32) | a;
}
#else
static uint64_t cmo_measure_start(void)
{
    uint64_t a, b;
    __asm__ __volatile__(
        "cpuid\n\t"
        "rdtsc\n\t"
        "mov %%eax, %0\n\t"
        "mov %%edx, %1" : "=r"(a), "=r"(b) : : "rax", "rbx", "rcx", "rdx");
    return ((uint64_t)b << 32) | a;
}

static uint64_t cmo_measure_end(void)
{
    uint32_t a, b;
    __asm__ __volatile__(
        "rdtscp\n\t"
        "mov %%eax, %0\n\t"
        "mov %%edx, %1\n\t"
        "cpuid" : "=r"(a), "=r"(b) : : "rax", "rbx", "rcx", "rdx");
    return ((uint64_t)b << 32) | a;
}
#endif

void cmo_init_prof(ProfileScratch* scr, int threadId, int entries)
{
    scr->threadId = threadId;
    scr->maxCount = entries;
    scr->entries = (ProfileEntry*)calloc(entries, sizeof(ProfileScratch));
    if (!scr->entries)
        assert(false && "Failed to allocaed ProfileEntrys");
}

void cmo_free_prof(ProfileScratch* scr)
{
    scr->maxCount = 0;
    free(scr->entries);
}

void cmo_prof_push_entry(ProfileScratch* scr, const char* name, 
                         uint64_t count, bool32_t exit)
{
    assert(scr->end < scr->maxCount - 1 && "Ran out of space on Profile stack");
    scr->entries[scr->end].func = name;
    scr->entries[scr->end].count = count;
    scr->entries[scr->end].exit = exit;
    ++scr->end;
}

void cmo_prof_start_timed_region(ProfileScratch* scr, const char* name)
{
    uint64_t count = cmo_measure_start();
    assert(scr->end < scr->maxCount - 1 && "Ran out of space on Profile stack");
    scr->entries[scr->end].func = name;
    scr->entries[scr->end].count = count;
    scr->entries[scr->end].exit = false;
    ++scr->end;
}

void cmo_prof_end_timed_region(ProfileScratch* scr, const char* name)
{
    uint64_t count = cmo_measure_end();
    assert(scr->end < scr->maxCount - 1 && "Ran out of space on Profile stack");
    scr->entries[scr->end].func = name;
    scr->entries[scr->end].count = count;
    scr->entries[scr->end].exit = true;
    ++scr->end;
}

void cmo_prof_print_file(ProfileScratch* scr, const char* fname)
{
    FILE* f = fopen(fname, "w");
    if (!f)
    {
        printf("Unable to open: %s\n", fname);
        return;
    }
    for (int i = 0; i < scr->end; ++i)
    {
        fprintf(f, "%s, %lu, %s", scr->entries[i].func, scr->entries[i].count, scr->entries[i].exit == true? "EXIT" : "ENTER");
        if (i != scr->end - 1)
        {
            fprintf(f, "\n");
        }
    }
    fflush(f);
    fclose(f);
}

void cmo_allocate_prof_array(int numThreads, int entries)
{
    printf("Allocating profiler space of %d threads\n", numThreads);
    profiler = calloc(numThreads, sizeof(ProfileScratch));
    for (int i = 0; i < numThreads; ++i)
    {
        cmo_init_prof(&profiler[i], i, entries);
    }
}

void cmo_free_prof_array(int numThreads)
{
    for (int i = 0; i < numThreads; ++i)
    {
        cmo_free_prof(&profiler[i]);
    }
    free(profiler);
}

// _Thread_local ProfileScratch threadProf;
_Thread_local int CmoProfileId;
// __thread int CmoProfileId;
// tss_t CmoProfileId;
// I think tss_t is borked everywhere due to lack of standards, what we can support as a backup is two options.
// Make sure that every threaded function has a threadId value, that can implicitly read by the macro
// OR
// re-engineer CmoProfileId to use pthread_key_t for a thread_local
// profiler/indexing array (the former may be harder to serialise).
// Can maybe just replace the thread_local with pthread_self
ProfileScratch* profiler;

#endif
#endif