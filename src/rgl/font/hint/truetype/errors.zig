//! What can go wrong running TrueType bytecode: FreeType's interpreter
//! errors (FreeType 2.14.3 `ttinterp.c`), with fontations' names where it
//! has one (`read-fonts/src/tables/glyf/bytecode/error.rs`, commit
//! 4f115482, MIT). Any of them ends the program, as in FreeType.

pub const Error = error{
    /// Operands past the end of the bytecode, or a call that ran off its
    /// end without ENDF (`Code_Overflow`).
    CodeOverflow,
    /// A jump into a program with no code (`Invalid_CodeRange`).
    InvalidCodeRange,
    /// An opcode with no instruction and no IDEF (`Invalid_Opcode`).
    UnhandledOpcode,
    /// A glyph program defined a function or instruction.
    DefinitionInGlyphProgram,
    /// FDEF or IDEF inside a definition.
    NestedDefinition,
    /// More definitions than `maxp` allows, or a bad number for one.
    TooManyDefinitions,
    /// ENDF with no call to return from.
    EndfWithoutCall,
    /// A call to a function never defined (`Invalid_Reference`).
    InvalidDefinition,
    /// Too few values on the stack, in pedantic mode.
    TooFewArguments,
    /// A push past the value stack, or a call past the call stack.
    StackOverflow,
    /// A bad index or zone, in pedantic mode (`Invalid_Reference`).
    InvalidReference,
    /// A bad jump, loop count, delta shift or the like (`Bad_Argument`).
    BadArgument,
    DivideByZero,
    /// The DEBUG instruction, which FreeType runs only in its debugger.
    DebugOpcode,
    /// Too many instructions, backward jumps or loop calls.
    ExecutionTooLong,
};
