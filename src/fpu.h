// x87 FPU state for RcppParallel workers.
//
// On Windows (x86/x64), R runs with the x87 FPU at 64-bit (extended)
// precision, but threads started by TBB begin at the Windows default of
// 53-bit precision. The mingw-w64 math library evaluates exp(), log() and
// related functions with x87 instructions, and long double arithmetic uses
// the x87 unit, so results depended on which thread computed them: the
// thread-invariance tests failed on Windows only.
//
// FpuState captures the control word of the thread that constructs a worker
// (R's main thread); FpuGuard applies it for the duration of a worker call
// and restores the previous one. Elsewhere both are empty no-ops.

#ifndef STS_FPU_H
#define STS_FPU_H

#if defined(_WIN32) && defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))

namespace sts_fpu {
inline unsigned short get_cw() {
  unsigned short cw;
  __asm__ __volatile__("fnstcw %0" : "=m"(cw));
  return cw;
}
inline void set_cw(unsigned short cw) {
  __asm__ __volatile__("fldcw %0" : : "m"(cw));
}
}  // namespace sts_fpu

struct FpuState {
  unsigned short cw;
  FpuState() : cw(sts_fpu::get_cw()) {}
};

class FpuGuard {
  unsigned short old_;
  bool changed_;
public:
  explicit FpuGuard(const FpuState& s)
    : old_(sts_fpu::get_cw()), changed_(old_ != s.cw) {
    if (changed_) sts_fpu::set_cw(s.cw);
  }
  ~FpuGuard() { if (changed_) sts_fpu::set_cw(old_); }
  FpuGuard(const FpuGuard&) = delete;
  FpuGuard& operator=(const FpuGuard&) = delete;
};

#else

struct FpuState {};
class FpuGuard {
public:
  explicit FpuGuard(const FpuState&) {}
};

#endif

#endif  // STS_FPU_H
