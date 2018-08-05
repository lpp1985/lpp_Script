from contextlib import contextmanager
import functools
import time
import warnings
import logging
import datetime

from pbsmrtpipe.constants import DEEP_DEBUG

log = logging.getLogger(__name__)


def timeit(func):

    @functools.wraps(func)
    def _(*args, **kwargs):
        started_at = datetime.datetime.now()
        r = func(*args, **kwargs)
        run_time = datetime.datetime.now() - started_at
        log.info("runtime for {f} {r} sec".format(r=run_time.seconds, f=func.__name__))
        return r

    return _


def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used."""

    def new_func(*args, **kwargs):
        warnings.warn("Call to deprecated function {}.".format(func.__name__),
                      category=DeprecationWarning)
        return func(*args, **kwargs)

    new_func.__name__ = func.__name__
    new_func.__doc__ = func.__doc__
    new_func.__dict__.update(func.__dict__)
    return new_func


def validate_type(*atype):
    """Simple validation of functions args

    Example:

    @validate_type(str)
    def my_func(name, message="Stuff"):
        print "Hello " + message + name

    """
    if len(atype) != 1:
        raise ValueError("Expected one arg. Got {n} args.".format(n=len(atype)))
    type_ = atype[0]

    def wrap(f):
        def wrapped_f(*args, **kw):
            if not isinstance(args[0], type_):
                raise TypeError("Expected type {t}. Got type {x} for {v}".format(t=type_, x=type(args[0]), v=args[0]))

            return f(*args)
        wrapped_f.__name = f.__name__
        wrapped_f.__doc__ = f.__doc__
        return wrapped_f
    return wrap


def validate_listable_type(*atype):
    """Validate a list of atype.

    @validate_listable_type(str)
        def example_func(a_list):
            return a_list

    @validate_listable_type(int)
    def example_int_func(a_list):
        return a_list

    """
    if len(atype) != 1:
        raise ValueError("Expected one arg. Got {n} args.".format(n=len(atype)))
    type_ = atype[0]

    def wrap(f):
        def wrapped_f(*args, **kw):
            for arg in args[0]:
                if not isinstance(arg, type_):
                    raise TypeError("Expected type {t}. Got type {x} for {v}.".format(t=type_, x=type(arg), v=args))

            return f(*args)
        return wrapped_f
    return wrap


@contextmanager
def ignored(*exceptions):
    try:
        yield
    except exceptions:
        pass


def trace(func):
    """
    Simple decorator to inspect the args, kwargs of a function.

    Useful debugging other decorators and functions.
    """
    if DEEP_DEBUG:
        log.debug("Called {f}".format(f=func.__name__))

        def wrap_func(*args, **kw):
            log.debug("args:")
            for a in args:
                log.debug(a)
            log.debug("kwargs:")
            for k, v in kw.iteritems():
                log.debug((k, v))
            result = func(*args, **kw)
            log.debug("Completed {f} with result {r}".format(f=func.__name__,
                                                             r=result))
            return result

        return wrap_func
    else:
        return func
