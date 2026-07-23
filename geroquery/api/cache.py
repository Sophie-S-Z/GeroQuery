"""Tiny version-aware in-process cache for expensive assemblies.

Keys include the data version, so a store rebuild transparently invalidates
everything. This is the MVP disk/Redis stand-in described in the PRD; the
interface (``get``/``set`` keyed by a tuple) would not change if swapped.
"""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Hashable
from typing import Any


class VersionedLRU:
    def __init__(self, maxsize: int = 512):
        self.maxsize = maxsize
        self._data: OrderedDict[Hashable, Any] = OrderedDict()

    def key(self, version: str, *parts: Hashable) -> Hashable:
        return (version, *parts)

    def get(self, key: Hashable):
        if key in self._data:
            self._data.move_to_end(key)
            return self._data[key]
        return None

    def set(self, key: Hashable, value: Any) -> None:
        self._data[key] = value
        self._data.move_to_end(key)
        while len(self._data) > self.maxsize:
            self._data.popitem(last=False)

    def clear(self) -> None:
        self._data.clear()
