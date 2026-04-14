#!/usr/bin/env python3
import os
import sys

project_root = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(project_root, "src")

if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

from main import main


if __name__ == "__main__":
    main()
