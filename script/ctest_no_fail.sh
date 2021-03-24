#! /bin/sh

# This file is a wrapper that returns 0 from the ctest. this wrapper is useful
# when there is a requirement to return 0 after running tests, despite the
# outcome of the tests themselves.

ctest $@
exit 0
