// Thin API client for the GeroQuery REST service.
const BASE = import.meta.env.VITE_API_BASE ?? "";

async function request(path, options) {
  const res = await fetch(`${BASE}${path}`, options);
  const body = await res.json().catch(() => ({}));
  if (!res.ok) {
    const msg = body?.error?.message ?? `Request failed (${res.status})`;
    throw new Error(msg);
  }
  return body;
}

const getJSON = (path) => request(path);

const postJSON = (path, payload) =>
  request(path, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });

/** The real NHANES clinical dataset. Its synthetic counterpart is
 *  `clinical_synthetic_csd` — deliberately a separate id so the UI can never
 *  render planted-effect numbers under a "real data" heading. */
export const REAL_CLINICAL_DATASET = "clinical_nhanes_slice";

export const api = {
  version: () => getJSON("/v1/version"),
  datasets: () => getJSON("/v1/datasets"),
  geneCard: (id, species) =>
    getJSON(`/v1/gene/${encodeURIComponent(id)}/card${species ? `?species=${species}` : ""}`),
  csd: (datasetId = REAL_CLINICAL_DATASET, nStrata = 6) =>
    postJSON("/v1/resilience/csd", { dataset_id: datasetId, n_strata: nStrata }),
};
