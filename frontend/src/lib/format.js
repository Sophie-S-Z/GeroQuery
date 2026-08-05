/** Number formatting. Every measurement is signed, fixed-width, and monospaced. */

import { format } from "d3-format";

const two = format("+.2f");
const three = format("+.3f");
const sci = format(".1e");

/** Signed effect size. The sign is information, so it is never dropped. */
export function fmtEffect(value, digits = 2) {
  if (value == null || Number.isNaN(value)) return "—";
  return digits === 3 ? three(value) : two(value);
}

/** A confidence interval, always shown next to the estimate it qualifies. */
export function fmtInterval(low, high) {
  if (low == null || high == null) return "—";
  return `[${two(low)}, ${two(high)}]`;
}

/**
 * p-values, with a floor rather than a fake zero.
 *
 * `0` is not a p-value any finite sample can produce, and printing it invites
 * a reader to treat a bound as a measurement.
 */
export function fmtP(value) {
  if (value == null || Number.isNaN(value)) return "—";
  if (value === 0) return "<1e-300";
  if (value < 0.001) return sci(value);
  return value.toFixed(3);
}

/** Hazard ratios are multiplicative, so they are unsigned and never "+1.44". */
export function fmtHR(value) {
  if (value == null || Number.isNaN(value)) return "—";
  return value.toFixed(3);
}

export function fmtHRInterval(low, high) {
  if (low == null || high == null) return "—";
  return `[${low.toFixed(3)}, ${high.toFixed(3)}]`;
}

export const VERDICT_LABEL = {
  increases: "rises with age",
  decreases: "falls with age",
  no_evidence: "no evidence",
};

/** The sentence a reader should take away, in the tool's own voice. */
export const VERDICT_SUB = {
  increases: "interval excludes zero",
  decreases: "interval excludes zero",
  no_evidence: "interval crosses zero",
};

export function fmtInt(value) {
  return value == null ? "—" : Number(value).toLocaleString();
}
