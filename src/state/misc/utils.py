"""
Misc utils
"""

import contextlib
import io


# -----------------------------------------------------------------------------
#   Functions
# -----------------------------------------------------------------------------


@contextlib.contextmanager
def buffered_stdout():
    """
    Use this to buffer print()'s until end of block.
    Useful when printing from processes, to keep entire print tnsrnm together.

    >>> with buffered_stdout():
    >>>     print("abc")
    >>>     print("123")
    """
    exc = None
    with io.StringIO() as buf, contextlib.redirect_stdout(buf):
        try:
            yield
        except Exception as e:
            # Catch exception, so curr output can be printed before raising the exception
            exc = e

        # noinspection PyUnresolvedReferences
        msg = buf.getvalue()

    print(msg, flush=True)

    # Raise the caught exception, if any
    if exc is not None:
        raise exc

    return
