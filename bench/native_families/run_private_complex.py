#!/usr/bin/env python3
"""Explicit complex-consumer qualification; reuses the private launch boundary."""
from complex_fixture_audit import ComplexFixture
from run_private_real import main

if __name__ == '__main__':
    raise SystemExit(main(ComplexFixture))
