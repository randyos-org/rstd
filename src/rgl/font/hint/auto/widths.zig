//! Stem widths and their sorting: FreeType 2.14.3 `src/autofit/aftypes.h`
//! (`AF_WidthRec`) and `afhints.c` (`af_sort_pos`,
//! `af_sort_and_quantize_widths`).

/// A width or position: in font units, scaled, and fitted to the grid.
pub const Width = struct {
    org: i32 = 0,
    cur: i32 = 0,
    fit: i32 = 0,
};

/// `af_sort_pos`: an insertion sort, stable.
pub fn sortPos(table: []i32) void {
    for (1..@max(table.len, 1)) |i| {
        var j = i;
        while (j > 0 and table[j] < table[j - 1]) : (j -= 1) {
            const t = table[j];
            table[j] = table[j - 1];
            table[j - 1] = t;
        }
    }
}

/// `af_sort_and_quantize_widths`: the first `n` of `table` sorted, each
/// cluster no wider than `threshold` replaced by a mean; the count left.
/// As FreeType's, a count of 0 becomes 1, keeping `table[0]`.
pub fn sortAndQuantize(table: []Width, n: usize, threshold: i32) usize {
    var count = n;
    if (count == 1) return 1;
    for (1..@max(count, 1)) |i| {
        var j = i;
        while (j > 0 and table[j].org < table[j - 1].org) : (j -= 1) {
            const t = table[j];
            table[j] = table[j - 1];
            table[j - 1] = t;
        }
    }
    var cur_idx: usize = 0;
    var cur_val = table[0].org;
    var i: usize = 1;
    while (i < count) : (i += 1) {
        if (table[i].org - cur_val > threshold or i == count - 1) {
            var sum: i32 = 0;
            // the last entry joins its cluster
            if (table[i].org - cur_val <= threshold and i == count - 1) i += 1;
            var j = cur_idx;
            while (j < i) : (j += 1) {
                sum += table[j].org;
                table[j].org = 0;
            }
            table[cur_idx].org = @divTrunc(sum, @as(i32, @intCast(j)));
            if (i < count - 1) {
                cur_idx = i + 1;
                cur_val = table[cur_idx].org;
            }
        }
    }
    var out: usize = 1;
    for (1..@max(count, 1)) |k| {
        if (table[k].org != 0) {
            table[out] = table[k];
            out += 1;
        }
    }
    count = out;
    return count;
}
