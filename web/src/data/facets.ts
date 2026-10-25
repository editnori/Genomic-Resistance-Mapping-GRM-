// Filter counts that predict: the count beside an option is the number of
// rows choosing it would return while every other control stays as it is.

/** One filter axis. `test` is null when the axis is at rest (no filter). */
export interface Axis<Row> {
  id: string;
  test: ((row: Row) => boolean) | null;
}

export const applyAxes = <Row>(rows: readonly Row[], axes: readonly Axis<Row>[]): Row[] =>
  rows.filter((row) => axes.every((a) => a.test === null || a.test(row)));

/** Rows that pass every axis except `lift`, the axis being counted. */
export const rowsWithout = <Row>(rows: readonly Row[], axes: readonly Axis<Row>[], lift: string): Row[] =>
  applyAxes(
    rows,
    axes.filter((a) => a.id !== lift),
  );

/** Count per value of a one-valued facet over the rows the other axes leave. */
export const countBy = <Row>(
  rows: readonly Row[],
  axes: readonly Axis<Row>[],
  lift: string,
  valueOf: (row: Row) => string,
): Map<string, number> => {
  const counts = new Map<string, number>();
  rowsWithout(rows, axes, lift).forEach((row) => {
    const v = valueOf(row);
    counts.set(v, (counts.get(v) ?? 0) + 1);
  });
  return counts;
};

/** Many-of selection: an empty set means every value. */
export const memberTest = <Row>(selected: ReadonlySet<string>, valueOf: (row: Row) => string): ((row: Row) => boolean) | null =>
  selected.size === 0 ? null : (row) => selected.has(valueOf(row));
