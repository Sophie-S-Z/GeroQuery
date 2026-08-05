/**
 * URL-as-state.
 *
 * The single most-requested thing missing from the Streamlit build was a
 * shareable link: you could not send someone a gene result. Here the query
 * string *is* the application state, so every view has an address, the back
 * button works, and a result can be cited.
 */

import { useCallback, useEffect, useState } from "react";

const DEFAULTS = { view: "gene", gene: "CDKN2A", species: "human", dir: "up" };

function read() {
  const params = new URLSearchParams(window.location.search);
  const state = { ...DEFAULTS };
  for (const key of Object.keys(DEFAULTS)) {
    const value = params.get(key);
    if (value) state[key] = value;
  }
  return state;
}

/**
 * @returns {[Record<string,string>, (patch: Record<string,string>) => void]}
 */
export function useUrlState() {
  const [state, setState] = useState(read);

  useEffect(() => {
    const onPop = () => setState(read());
    window.addEventListener("popstate", onPop);
    return () => window.removeEventListener("popstate", onPop);
  }, []);

  const update = useCallback((patch) => {
    setState((previous) => {
      const next = { ...previous, ...patch };
      const params = new URLSearchParams();
      for (const [key, value] of Object.entries(next)) {
        // Omit anything at its default so a shared link carries only what the
        // sender actually chose. A URL full of redundant defaults is noise.
        if (value && value !== DEFAULTS[key]) params.set(key, value);
      }
      const search = params.toString();
      window.history.pushState(
        null,
        "",
        search ? `${window.location.pathname}?${search}` : window.location.pathname,
      );
      return next;
    });
  }, []);

  return [state, update];
}
