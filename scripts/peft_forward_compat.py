"""Adapt PEFT's optional keyword arguments to a custom backbone signature."""
import inspect


def adapt_forward_for_peft(model):
    model_class = type(model)
    original = model_class.forward  # unbound; dispatch to the current replica
    parameters = inspect.signature(original).parameters
    if any(p.kind == inspect.Parameter.VAR_KEYWORD for p in parameters.values()):
        return model

    def compatible_forward(self, *args, **kwargs):
        forwarded = {}
        for name, value in kwargs.items():
            if name in parameters:
                forwarded[name] = value
            elif value is not None:
                raise TypeError(
                    f"Backbone does not accept {name!r} with a non-None value; "
                    "this adapter check uses unpadded, fixed-length sequences."
                )
        return original(self, *args, **forwarded)

    model_class.forward = compatible_forward
    return model
