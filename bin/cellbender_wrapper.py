#!/usr/bin/env python3
"""
cellbender_wrapper.py — FIX-22b
Patches CellBender 0.3.0 for PyTorch 2.4+ checkpoint compatibility.

Problem:
  PyTorch 2.4 changed nn.Module internals to use weakref-backed containers
  in hooks, Pyro guide internals, and other deep model attributes.
  CellBender's checkpoint.py calls torch.save(model_obj, ...) which tries
  to pickle the entire model graph.  This fails with:
      TypeError: cannot pickle 'weakref.ReferenceType' object

  Without the checkpoint tarball, the posterior computation step crashes:
      AssertionError: Checkpoint file ckpt.tar.gz does not exist

  FIX-22 attempted to strip hook dicts, but the weakrefs live deeper in
  Pyro's variational inference internals, not just in PyTorch hooks.

Fix (FIX-22b):
  Use copyreg to register a global pickle handler for weakref.ref objects.
  When pickle encounters ANY weakref anywhere in the object graph, it
  serialises it as None instead of crashing.  This is safe because:
    1. The checkpoint is saved AFTER training completes
    2. Hooks and weakref-backed observers are only needed during training
    3. The loaded checkpoint is only used for posterior inference (forward
       passes), which doesn't require hook machinery
  Additionally, still strip hook dicts as a belt-and-suspenders measure.

Usage (drop-in replacement for the `cellbender` CLI):
  python3 cellbender_wrapper.py remove-background --input ... --output ...

Place in bin/ alongside the other pipeline scripts.
"""

import sys
import copyreg
import weakref
import torch
import torch.nn as nn

# ---------------------------------------------------------------------------
# FIX-22b: Make weakref.ref globally picklable
# ---------------------------------------------------------------------------
# When pickle encounters a weakref.ref object anywhere in the object graph,
# it will call _pickle_weakref to serialize it.  On deserialization,
# _restore_as_none returns None (the weakref target is no longer needed).

def _restore_as_none():
    """Unpickle placeholder: returns None where a weakref used to be."""
    return None


def _pickle_weakref(wr):
    """Reduce a weakref.ref to a picklable form (restores as None)."""
    return (_restore_as_none, ())


# Register for weakref.ref (covers weakref.ReferenceType)
copyreg.pickle(weakref.ref, _pickle_weakref)

# Also handle weakref.CallableProxyType if encountered
try:
    copyreg.pickle(weakref.CallableProxyType, _pickle_weakref)
except (TypeError, AttributeError):
    pass  # Not all Python versions expose this as a registerable type

print("FIX-22b: copyreg weakref.ref pickle handler registered", flush=True)


# ---------------------------------------------------------------------------
# Belt-and-suspenders: also monkey-patch torch.save to strip hooks
# ---------------------------------------------------------------------------
_original_torch_save = torch.save

_HOOK_ATTRS = (
    "_backward_hooks",
    "_forward_hooks",
    "_forward_pre_hooks",
    "_backward_pre_hooks",
    "_state_dict_hooks",
    "_load_state_dict_pre_hooks",
    "_state_dict_pre_hooks",
    "_load_state_dict_post_hooks",
    # PyTorch 2.4+ additions
    "_forward_hooks_with_kwargs",
    "_forward_hooks_always_called",
    "_forward_pre_hooks_with_kwargs",
)


def _strip_hooks(module: nn.Module):
    """Replace hook containers with plain dicts to remove weakref wrappers."""
    for m in module.modules():
        for attr in _HOOK_ATTRS:
            if hasattr(m, attr):
                try:
                    val = getattr(m, attr)
                    if isinstance(val, dict):
                        setattr(m, attr, dict(val))
                    else:
                        setattr(m, attr, {})
                except Exception:
                    try:
                        setattr(m, attr, {})
                    except Exception:
                        pass


def _safe_torch_save(obj, f, *args, **kwargs):
    """torch.save wrapper: strips hooks from nn.Modules before saving.

    Combined with the copyreg weakref handler, this handles both
    hook-level and deep-model weakrefs.
    """
    if isinstance(obj, nn.Module):
        _strip_hooks(obj)

    try:
        _original_torch_save(obj, f, *args, **kwargs)
    except TypeError as exc:
        if "weakref" not in str(exc) and "cannot pickle" not in str(exc):
            raise
        # copyreg should have handled this, but as a last resort
        # try saving just the state_dict
        print(f"FIX-22b: Full model save failed ({exc}), attempting state_dict fallback...",
              flush=True)
        if isinstance(obj, nn.Module):
            _original_torch_save(obj.state_dict(), f, *args, **kwargs)
            print("FIX-22b: state_dict fallback save succeeded", flush=True)
        else:
            raise


# Apply the patch
torch.save = _safe_torch_save
print("FIX-22b: torch.save patched (hook strip + copyreg weakref handler)", flush=True)

# ---------------------------------------------------------------------------
# Run CellBender
# ---------------------------------------------------------------------------
from cellbender.base_cli import main  # noqa: E402

sys.exit(main())
