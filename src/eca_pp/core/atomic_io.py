"""Atomic file writes — kill the 0-byte / half-written-file failure class at the source.

Writes go to a temp file in the SAME directory, then ``os.replace`` renames it into
place (atomic on POSIX within one filesystem): a reader sees either the old complete
file or the new one, never a torn one. A crash leaves only the ``.tmp``, which no
output glob can mistake for a finished file.

v0.1 carries only the pieces the standardize CLI needs (result.json); the h5ad and
figure writers join with F7/F8.
"""

import os
from contextlib import contextmanager


@contextmanager
def atomic_write(path):
    """Yield a temp path next to ``path``; atomically rename onto ``path`` on success,
    remove the temp and re-raise on error. Creates parent dirs."""
    path = os.fspath(path)
    dirname = os.path.dirname(path) or "."
    os.makedirs(dirname, exist_ok=True)
    tmp = f"{path}.{os.getpid()}.tmp"
    try:
        yield tmp
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass  # temp may not exist yet; never mask the real error
        raise


def write_bytes_atomic(path, data: bytes) -> None:
    """Write ``data`` to ``path`` atomically."""
    with atomic_write(path) as tmp:
        with open(tmp, "wb") as fh:
            fh.write(data)
