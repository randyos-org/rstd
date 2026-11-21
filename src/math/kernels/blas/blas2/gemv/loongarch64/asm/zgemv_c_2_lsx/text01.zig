//! Section 1 of 2 of the text of `zgemv_c`, which `__root__.zig` joins, from OpenBLAS `kernel/loongarch64/zgemv_t_2_lsx.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .altmacro
    \\ .macro push_if_used regs, fregs
    \\ .if \regs > 0
    \\ addi.d $sp, $sp, -(\regs << 3)
    \\ push_regs 0, \regs - 1
    \\ .endif
    \\ .if \fregs > 0
    \\ addi.d $sp, $sp, -(\fregs << 3)
    \\ push_fregs 0, \fregs - 1
    \\ .endif
    \\ .endm
    \\ .macro pop_if_used regs, fregs
    \\ .if \fregs > 0
    \\ pop_fregs 0, \fregs - 1
    \\ addi.d $sp, $sp, \fregs << 3
    \\ .endif
    \\ .if \regs > 0
    \\ pop_regs 0, \regs - 1
    \\ addi.d $sp, $sp, \regs << 3
    \\ .endif
    \\ .endm
    \\ .macro push_regs from, to
    \\ .if \to >= 0
    \\ st.d $s0, $sp, 0 << 3
    \\ .endif
    \\ .if \to >= 1
    \\ st.d $s1, $sp, 1 << 3
    \\ .endif
    \\ .if \to >= 2
    \\ st.d $s2, $sp, 2 << 3
    \\ .endif
    \\ .if \to >= 3
    \\ st.d $s3, $sp, 3 << 3
    \\ .endif
    \\ .if \to >= 4
    \\ st.d $s4, $sp, 4 << 3
    \\ .endif
    \\ .if \to >= 5
    \\ st.d $s5, $sp, 5 << 3
    \\ .endif
    \\ .if \to >= 6
    \\ st.d $s6, $sp, 6 << 3
    \\ .endif
    \\ .if \to >= 7
    \\ st.d $s7, $sp, 7 << 3
    \\ .endif
    \\ .if \to >= 8
    \\ st.d $s8, $sp, 8 << 3
    \\ .endif
    \\ .endm
    \\ .macro pop_regs from, to
    \\ .if \to >= 0
    \\ ld.d $s0, $sp, 0 << 3
    \\ .endif
    \\ .if \to >= 1
    \\ ld.d $s1, $sp, 1 << 3
    \\ .endif
    \\ .if \to >= 2
    \\ ld.d $s2, $sp, 2 << 3
    \\ .endif
    \\ .if \to >= 3
    \\ ld.d $s3, $sp, 3 << 3
    \\ .endif
    \\ .if \to >= 4
    \\ ld.d $s4, $sp, 4 << 3
    \\ .endif
    \\ .if \to >= 5
    \\ ld.d $s5, $sp, 5 << 3
    \\ .endif
    \\ .if \to >= 6
    \\ ld.d $s6, $sp, 6 << 3
    \\ .endif
    \\ .if \to >= 7
    \\ ld.d $s7, $sp, 7 << 3
    \\ .endif
    \\ .if \to >= 8
    \\ ld.d $s8, $sp, 8 << 3
    \\ .endif
    \\ .endm
    \\ .macro push_fregs from, to
    \\ .if \to >= 0
    \\ fst.d $fs0, $sp, 0 << 3
    \\ .endif
    \\ .if \to >= 1
    \\ fst.d $fs1, $sp, 1 << 3
    \\ .endif
    \\ .if \to >= 2
    \\ fst.d $fs2, $sp, 2 << 3
    \\ .endif
    \\ .if \to >= 3
    \\ fst.d $fs3, $sp, 3 << 3
    \\ .endif
    \\ .if \to >= 4
    \\ fst.d $fs4, $sp, 4 << 3
    \\ .endif
    \\ .if \to >= 5
    \\ fst.d $fs5, $sp, 5 << 3
    \\ .endif
    \\ .if \to >= 6
    \\ fst.d $fs6, $sp, 6 << 3
    \\ .endif
    \\ .if \to >= 7
    \\ fst.d $fs7, $sp, 7 << 3
    \\ .endif
    \\ .endm
    \\ .macro pop_fregs from, to
    \\ .if \to >= 0
    \\ fld.d $fs0, $sp, 0 << 3
    \\ .endif
    \\ .if \to >= 1
    \\ fld.d $fs1, $sp, 1 << 3
    \\ .endif
    \\ .if \to >= 2
    \\ fld.d $fs2, $sp, 2 << 3
    \\ .endif
    \\ .if \to >= 3
    \\ fld.d $fs3, $sp, 3 << 3
    \\ .endif
    \\ .if \to >= 4
    \\ fld.d $fs4, $sp, 4 << 3
    \\ .endif
    \\ .if \to >= 5
    \\ fld.d $fs5, $sp, 5 << 3
    \\ .endif
    \\ .if \to >= 6
    \\ fld.d $fs6, $sp, 6 << 3
    \\ .endif
    \\ .if \to >= 7
    \\ fld.d $fs7, $sp, 7 << 3
    \\ .endif
    \\ .endm
    \\ .macro GLD pre_op:req, suf_op=0, out:req, src:req, offset:req , more:vararg
    \\ .ifeqs "\suf_op", "0"
    \\ \pre_op\()ld \out, \src, \offset
    \\ .else
    \\ \pre_op\()ld.\suf_op \out, \src, \offset
    \\ .endif
    \\ .ifnb \more
    \\ GLD \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro GLD_INC pre_op:req, suf_op=0, inc:req, out:req, src:req, offset:req , more:vararg
    \\ .ifeqs "\suf_op", "0"
    \\ \pre_op\()ld \out, \src, \offset
    \\ .else
    \\ \pre_op\()ld.\suf_op \out, \src, \offset
    \\ .endif
    \\ addi.d \src, \src, \inc
    \\ .ifnb \more
    \\ GLD_INC \pre_op, \suf_op, \inc, \more
    \\ .endif
    \\ .endm
    \\ .macro GLDX pre_op:req, suf_op=0, out:req, src:req, offset:req , more:vararg
    \\ .ifeqs "\suf_op", "0"
    \\ \pre_op\()ldx \out, \src, \offset
    \\ .else
    \\ \pre_op\()ldx.\suf_op \out, \src, \offset
    \\ .endif
    \\ .ifnb \more
    \\ GLDX \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro GLDREPL pre_op:req, suf_op:req, out:req, src:req, offset:req , more:vararg
    \\ \pre_op\()ldrepl.\suf_op \out, \src, \offset
    \\ .ifnb \more
    \\ GLDREPL \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro GST pre_op:req, suf_op=0, src:req, dst:req, offset:req , more:vararg
    \\ .ifeqs "\suf_op", "0"
    \\ \pre_op\()st \src, \dst, \offset
    \\ .else
    \\ \pre_op\()st.\suf_op \src, \dst, \offset
    \\ .endif
    \\ .ifnb \more
    \\ GST \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro GMUL pre_op, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \\ \pre_op\()mul.\suf_op \out, \in0, \in1
    \\ .ifnb \more
    \\ GMUL \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro GMADD pre_op, suf_op:req, out:req, in0:req, in1:req, in2:req, more:vararg
    \\ \pre_op\()madd.\suf_op \out, \in0, \in1, \in2
    \\ .ifnb \more
    \\ GMADD \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro GADD pre_op, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \\ \pre_op\()add.\suf_op \out, \in0, \in1
    \\ .ifnb \more
    \\ GADD \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro GADDI pre_op, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \\ \pre_op\()addi.\suf_op \out, \in0, \in1
    \\ .ifnb \more
    \\ GADDI \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro GSUB pre_op, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \\ \pre_op\()sub.\suf_op \out, \in0, \in1
    \\ .ifnb \more
    \\ GSUB \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro GSLLI pre_op, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \\ \pre_op\()slli.\suf_op \out, \in0, \in1
    \\ .ifnb \more
    \\ GSLLI \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro GINSVE0 pre_op:req, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \\ \pre_op\()insve0.\suf_op \out, \in0, \in1
    \\ .ifnb \more
    \\ GINSVE0 \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro GXOR pre_op:req, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \\ .ifnb \pre_op
    \\ \pre_op\()xor.v \out, \in0, \in1
    \\ .else
    \\ xor.\suf_op \out, \in0, \in1
    \\ .endif
    \\ .ifnb \more
    \\ GXOR \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro GPERMI pre_op:req, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \\ \pre_op\()permi.\suf_op \out, \in0, \in1
    \\ .ifnb \more
    \\ GPERMI \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro GNMSUB pre_op:req, suf_op:req, out:req, in0:req, in1:req, in2:req, more:vararg
    \\ \pre_op\()nmsub.\suf_op \out, \in0, \in1, \in2
    \\ .ifnb \more
    \\ GNMSUB \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro GPRELD in0:req, in1:req, in2:req, more:vararg
    \\ preld \in0, \in1, \in2
    \\ .ifnb \more
    \\ GPRELD \more
    \\ .endif
    \\ .endm
    \\ .macro GPACKEV pre_op:req, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \\ \pre_op\()packev.\suf_op \out, \in0, \in1
    \\ .ifnb \more
    \\ GPACKEV \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro GPACKOD pre_op:req, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \\ \pre_op\()packod.\suf_op \out, \in0, \in1
    \\ .ifnb \more
    \\ GPACKOD \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro GSHUF4I pre_op:req, suf_op:req, out:req, in0:req, in1:req , more:vararg
    \\ \pre_op\()shuf4i.\suf_op \out, \in0, \in1
    \\ .ifnb \more
    \\ GSHUF4I \pre_op, \suf_op, \more
    \\ .endif
    \\ .endm
    \\ .macro TRANSF2G name, pre_op:req, suf_op:req, more:vararg
    \\ .ifeqs "\pre_op\()\suf_op", "vfs"
    \\ \name v, w, \more
    \\ .endif
    \\ .ifeqs "\pre_op\()\suf_op", "vfd"
    \\ \name v, d, \more
    \\ .endif
    \\ .ifeqs "\pre_op\()\suf_op", "xvfs"
    \\ \name xv, w, \more
    \\ .endif
    \\ .ifeqs "\pre_op\()\suf_op", "xvfd"
    \\ \name xv, d, \more
    \\ .endif
    \\ .endm
    \\
;
