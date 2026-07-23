// Thin API client for the GeroQuery REST service.
const BASE = import.meta.env.VITE_API_BASE ?? "";

async function getJSON(path) {
  const res = await fetch(`${BASE}${path}`);
  const body = await res.json().catch(() => ({}));
  if (!res.ok) {
    const msg = body?.error?.message ?? `Request failed (${res.status})`;
    throw new Error(msg);
  }
  return body;
}

export const api = {
  version: () => getJSON("/v1/version"),
  geneCard: (id, species) =>
    getJSON(`/v1/gene/${encodeURIComponent(id)}/card${species ? `?species=${species}` : ""}`),
};
