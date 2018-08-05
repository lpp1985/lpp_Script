
import platform

def add_override_distribute_option(p):
    g = p.add_mutually_exclusive_group()
    g.add_argument('--force-distributed', action='store_const', const=True, default=None,
                   help="Override XML settings to enable distributed mode (if cluster manager is provided)")
    g.add_argument('--local-only', action='store_const', const=True, default=None,
                   help="Override XML settings to disable distributed mode. All Task will be submitted to {n}".format(n=platform.node()))
    return p


def add_override_chunked_mode(p):
    g = p.add_mutually_exclusive_group()
    g.add_argument('--force-chunk-mode', action='store_const', const=True, default=None, help="Override to enable Chunk mode")
    g.add_argument('--disable-chunk-mode', action='store_const', const=True, default=None, help="Override to disable Chunk mode")
    return p
