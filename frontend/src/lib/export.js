/**
 * Rows out of the page as CSV, for a spreadsheet or a supplement.
 *
 * Writes exactly the values the page shows; nothing is recomputed on the way
 * out, so a downloaded table and the page cannot disagree. Every cell comes
 * from the shipped data, never from pasted text, so no cell can carry a
 * spreadsheet formula a reader did not put there.
 */

function cell(value) {
  if (value == null) return "";
  const text = String(value);
  return /[",\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

/** Download `rows` as CSV, with `columns` naming the fields in order. */
export function downloadCsv(filename, columns, rows) {
  const lines = [columns.join(","), ...rows.map((row) => columns.map((c) => cell(row[c])).join(","))];
  const url = URL.createObjectURL(new Blob([`${lines.join("\n")}\n`], { type: "text/csv" }));
  const link = Object.assign(document.createElement("a"), { href: url, download: filename });
  link.click();
  URL.revokeObjectURL(url);
}
