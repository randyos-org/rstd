//! One instruction: FreeType 2.14.3 `ttinterp.c` `TT_RunIns`'s switch.
//! The instructions that move points (#807) are not here yet and report
//! `error.UnhandledOpcode`.

const opcode = @import("opcode.zig");
const Engine = @import("Engine.zig");
const stack = @import("ops/stack.zig");
const arith = @import("ops/arith.zig");
const flow = @import("ops/flow.zig");
const calls = @import("ops/calls.zig");
const state = @import("ops/state.zig");
const data = @import("ops/data.zig");
const tags = @import("ops/tags.zig");

pub fn step(e: *Engine) Engine.Error!void {
    const op = e.op;
    if (op >= opcode.pushb_first and op < opcode.mdrp_first) return stack.push(e);
    // MDRP and MIRP move points (#807)
    if (op >= opcode.mdrp_first) return error.UnhandledOpcode;
    switch (@as(opcode.Opcode, @fromBackingInt(op))) {
        .svtca0, .svtca1, .spvtca0, .spvtca1, .sfvtca0, .sfvtca1 => state.toAxis(e),
        .spvtl0, .spvtl1, .sfvtl0, .sfvtl1 => try state.toLine(e),
        .spvfs, .sfvfs => state.fromStack(e),
        .gpv, .gfv => state.get(e),
        .sfvtpv => state.freedomToProjection(e),
        .sdpvtl0, .sdpvtl1 => try state.dualToLine(e),
        .srp0 => state.referencePoint(e, 0),
        .srp1 => state.referencePoint(e, 1),
        .srp2 => state.referencePoint(e, 2),
        .szp0 => try state.zonePointer(e, 0),
        .szp1 => try state.zonePointer(e, 1),
        .szp2 => try state.zonePointer(e, 2),
        .szps => try state.zonePointer(e, 3),
        .sloop => try state.setLoop(e),
        .rtg => state.setRound(e, .grid),
        .rthg => state.setRound(e, .half_grid),
        .rtdg => state.setRound(e, .double_grid),
        .rdtg => state.setRound(e, .down_to_grid),
        .rutg => state.setRound(e, .up_to_grid),
        .roff => state.setRound(e, .off),
        .sround => state.superRound(e, .super),
        .s45round => state.superRound(e, .super45),
        .smd, .scvtci, .sswci => state.setDistance(e),
        .ssw => state.singleWidth(e),
        .flipon => e.gs.auto_flip = true,
        .flipoff => e.gs.auto_flip = false,
        .sdb => state.deltaBase(e),
        .sds => try state.deltaShift(e),
        .instctrl => try state.instructControl(e),
        .scanctrl => state.scanControl(e),
        .scantype => state.scanType(e),
        // SANGW and AA are no longer supported; they pop their argument
        .sangw, .aa, .pop, .eif, .nround00, .nround01, .nround10, .nround11 => {},
        .debug => return error.DebugOpcode,
        .@"else" => try flow.@"else"(e),
        .@"if" => try flow.@"if"(e),
        .jmpr => try flow.jump(e),
        .jrot => try flow.jumpIf(e, true),
        .jrof => try flow.jumpIf(e, false),
        .dup => stack.dup(e),
        .clear => stack.clear(e),
        .swap => stack.swap(e),
        .depth => stack.depth(e),
        .cindex => try stack.cindex(e),
        .mindex => try stack.mindex(e),
        .roll => stack.roll(e),
        .npushb, .npushw => try stack.push(e),
        .loopcall => try calls.loopcall(e),
        .call => try calls.call(e),
        .fdef => try calls.define(e, false),
        .idef => try calls.define(e, true),
        .endf => try calls.endf(e),
        .ws => try data.writeStorage(e),
        .rs => try data.readStorage(e),
        .wcvtp => try data.writeCvtPixels(e),
        .wcvtf => try data.writeCvtFunits(e),
        .rcvt => try data.readCvt(e),
        .deltac1, .deltac2, .deltac3 => try data.deltaCvt(e),
        .gc0, .gc1 => try data.getCoordinate(e),
        .md0, .md1 => try data.measureDistance(e),
        .mppem => data.ppem(e),
        .mps => data.pointSize(e),
        .getinfo => data.getInfo(e),
        .getvariation => try data.getVariation(e),
        .getdata => try data.getData(e),
        .lt => try arith.binary(e, .lt),
        .lteq => try arith.binary(e, .lteq),
        .gt => try arith.binary(e, .gt),
        .gteq => try arith.binary(e, .gteq),
        .eq => try arith.binary(e, .eq),
        .neq => try arith.binary(e, .neq),
        .@"and" => try arith.binary(e, .@"and"),
        .@"or" => try arith.binary(e, .@"or"),
        .add => try arith.binary(e, .add),
        .sub => try arith.binary(e, .sub),
        .div => try arith.binary(e, .div),
        .mul => try arith.binary(e, .mul),
        .max => try arith.binary(e, .max),
        .min => try arith.binary(e, .min),
        .abs => arith.unary(e, .abs),
        .neg => arith.unary(e, .neg),
        .floor => arith.unary(e, .floor),
        .ceiling => arith.unary(e, .ceiling),
        .not => arith.unary(e, .not),
        .odd => arith.unary(e, .odd),
        .even => arith.unary(e, .even),
        .round00, .round01, .round10, .round11 => arith.unary(e, .round),
        .flippt => try tags.flipPoints(e),
        .fliprgon => try tags.flipRange(e, true),
        .fliprgoff => try tags.flipRange(e, false),
        .utp => try tags.untouch(e),
        // the instructions that move points (#807)
        .isect, .alignpts, .mdap0, .mdap1, .iup0, .iup1, .shp0, .shp1, .shc0, .shc1, .shz0, .shz1, .shpix, .ip, .msirp0, .msirp1, .alignrp, .miap0, .miap1, .scfs, .deltap1, .deltap2, .deltap3 => return error.UnhandledOpcode,
        _ => try calls.unknown(e),
    }
}
