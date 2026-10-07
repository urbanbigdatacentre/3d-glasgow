"""
Minimal ctypes binding to the meshoptimizer C library (vendor/libmeshoptimizer.dylib),
covering what the tile writer needs for EXT_meshopt_compression:
vertex cache / fetch optimisation, vertex + index buffer encoding.

Build the library with tools/build_meshopt.sh (clang++, ~10 s).
"""
import ctypes
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_LIBNAME = "libmeshoptimizer.dylib" if sys.platform == "darwin" else "libmeshoptimizer.so"
_LIBPATH = os.path.join(_HERE, "vendor", _LIBNAME)
if not os.path.exists(_LIBPATH):
    raise ImportError("meshoptimizer library not found at %s; run tools/build_meshopt.sh" % _LIBPATH)
_lib = ctypes.CDLL(_LIBPATH)

c_size_t, c_void_p, c_int = ctypes.c_size_t, ctypes.c_void_p, ctypes.c_int
_lib.meshopt_encodeVertexVersion.argtypes = [c_int]
_lib.meshopt_encodeIndexVersion.argtypes = [c_int]
_lib.meshopt_encodeVertexBufferBound.restype = c_size_t
_lib.meshopt_encodeVertexBufferBound.argtypes = [c_size_t, c_size_t]
_lib.meshopt_encodeVertexBuffer.restype = c_size_t
_lib.meshopt_encodeVertexBuffer.argtypes = [c_void_p, c_size_t, c_void_p, c_size_t, c_size_t]
_lib.meshopt_encodeIndexBufferBound.restype = c_size_t
_lib.meshopt_encodeIndexBufferBound.argtypes = [c_size_t, c_size_t]
_lib.meshopt_encodeIndexBuffer.restype = c_size_t
_lib.meshopt_encodeIndexBuffer.argtypes = [c_void_p, c_size_t, c_void_p, c_size_t]
_lib.meshopt_optimizeVertexCache.argtypes = [c_void_p, c_void_p, c_size_t, c_size_t]
_lib.meshopt_optimizeVertexFetch.restype = c_size_t
_lib.meshopt_optimizeVertexFetch.argtypes = [c_void_p, c_void_p, c_size_t, c_void_p, c_size_t, c_size_t]

_lib.meshopt_simplify.restype = c_size_t
_lib.meshopt_simplify.argtypes = [c_void_p, c_void_p, c_size_t, c_void_p, c_size_t, c_size_t,
                                  c_size_t, ctypes.c_float, ctypes.c_uint, ctypes.POINTER(ctypes.c_float)]
SIMPLIFY_LOCK_BORDER = 1 << 0
SIMPLIFY_ERROR_ABSOLUTE = 1 << 2


def simplify(positions, indices, target_error, target_index_count=0, lock_border=False):
    """
    Error-bounded simplification. positions: float32 (N,3) in metres; indices: uint32 (M,).
    target_error is an absolute distance in metres (meshopt_SimplifyErrorAbsolute).
    Returns (new_indices uint32, result_error_m). Vertices are not compacted.
    """
    positions = np.ascontiguousarray(positions, dtype=np.float32)
    indices = np.ascontiguousarray(indices, dtype=np.uint32)
    dst = np.empty_like(indices)
    err = ctypes.c_float(0.0)
    opts = SIMPLIFY_ERROR_ABSOLUTE | (SIMPLIFY_LOCK_BORDER if lock_border else 0)
    n = _lib.meshopt_simplify(_ptr(dst), _ptr(indices), indices.size, _ptr(positions), len(positions), 12,
                              target_index_count, ctypes.c_float(target_error), opts, ctypes.byref(err))
    return dst[:n], float(err.value)


# Vertex codec v0 for the widest decoder compatibility (CesiumJS bundles the JS decoder);
# index codec v1 has been the default since 2020.
_lib.meshopt_encodeVertexVersion(0)
_lib.meshopt_encodeIndexVersion(1)


def _ptr(a):
    return a.ctypes.data_as(c_void_p)


def optimize(vertices, indices):
    """
    vertices: uint8 (N, stride) C-contiguous interleaved vertex records
    indices : uint32 (M,) triangle list
    Returns (vertices_reordered, indices_remapped) after vertex-cache and vertex-fetch
    optimisation, which is what makes the encoders effective.
    """
    vertices = np.ascontiguousarray(vertices, dtype=np.uint8)
    indices = np.ascontiguousarray(indices, dtype=np.uint32)
    n_vert, stride = vertices.shape
    idx_opt = np.empty_like(indices)
    _lib.meshopt_optimizeVertexCache(_ptr(idx_opt), _ptr(indices), indices.size, n_vert)
    vert_opt = np.empty_like(vertices)
    unique = _lib.meshopt_optimizeVertexFetch(_ptr(vert_opt), _ptr(idx_opt), idx_opt.size, _ptr(vertices), n_vert, stride)
    return vert_opt[:unique], idx_opt


def encode_vertex_buffer(vertices):
    """vertices: uint8 (N, stride). Returns bytes (EXT_meshopt_compression mode ATTRIBUTES)."""
    vertices = np.ascontiguousarray(vertices, dtype=np.uint8)
    n_vert, stride = vertices.shape
    bound = _lib.meshopt_encodeVertexBufferBound(n_vert, stride)
    buf = np.empty(bound, dtype=np.uint8)
    n = _lib.meshopt_encodeVertexBuffer(_ptr(buf), bound, _ptr(vertices), n_vert, stride)
    if n == 0:
        raise RuntimeError("meshopt_encodeVertexBuffer failed")
    return buf[:n].tobytes()


def encode_index_buffer(indices, n_vert):
    """indices: uint32 (M,) triangle list. Returns bytes (mode TRIANGLES)."""
    indices = np.ascontiguousarray(indices, dtype=np.uint32)
    bound = _lib.meshopt_encodeIndexBufferBound(indices.size, n_vert)
    buf = np.empty(bound, dtype=np.uint8)
    n = _lib.meshopt_encodeIndexBuffer(_ptr(buf), bound, _ptr(indices), indices.size)
    if n == 0:
        raise RuntimeError("meshopt_encodeIndexBuffer failed")
    return buf[:n].tobytes()


if __name__ == "__main__":
    # smoke test on a grid mesh
    rng = np.random.default_rng(0)
    g = 200
    xs, ys = np.meshgrid(np.arange(g), np.arange(g))
    pos = np.stack([xs.ravel() * 10, (rng.random(g * g) * 300).astype(int), ys.ravel() * 10], axis=1).astype(np.int16)
    fid = (np.arange(g * g) // 97).astype(np.uint16)
    rec = np.zeros((g * g, 12), dtype=np.uint8)
    rec[:, 0:6] = pos.view(np.uint8).reshape(-1, 6)
    rec[:, 8:10] = fid.view(np.uint8).reshape(-1, 2)
    i = (np.arange(g - 1)[:, None] * g + np.arange(g - 1)[None, :]).ravel()
    tris = np.concatenate([np.stack([i, i + 1, i + g], 1), np.stack([i + 1, i + g + 1, i + g], 1)]).astype(np.uint32).ravel()
    v2, i2 = optimize(rec, tris)
    ev = encode_vertex_buffer(v2)
    ei = encode_index_buffer(i2, len(v2))
    print("vertices %d x %dB = %d B -> %d B (%.2f B/vertex)" % (len(v2), 12, len(v2) * 12, len(ev), len(ev) / len(v2)))
    print("indices  %d (%d tris) = %d B -> %d B (%.2f B/tri)" % (i2.size, i2.size // 3, i2.size * 2, len(ei), len(ei) / (i2.size // 3)))
    assert sorted(np.ascontiguousarray(v2[:, 8:10]).view(np.uint16).ravel().tolist()) == sorted(fid.tolist())
    print("ok")
