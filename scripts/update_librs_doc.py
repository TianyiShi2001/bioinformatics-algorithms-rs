# Copyright (C) 2020 Tianyi Shi
#
# This file is part of rust-bio-edu.
#
# rust-bio-edu is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# rust-bio-edu is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with rust-bio-edu.  If not, see <http://www.gnu.org/licenses/>.

"""Update the content of README.md to src/lib.rs
"""

import os

os.chdir(os.path.dirname(os.path.dirname(__file__)))

src = "README.md"
dst = "src/lib.rs"

if __name__ == "__main__":
    readme_formatted = ""
    with open(src) as f:
        for line in f.readlines():
            readme_formatted += "//! "
            readme_formatted += line
    res = ""
    with open(dst) as f:
        before = True
        deleting = False
        inserted = False
        lines = f.readlines()
        i = 0
        while lines[i][:3] != "//!":
            res += lines[i]
            i += 1
        while lines[i][:3] == "//!":
            i += 1
        res += readme_formatted
        while i < len(lines):
            res += lines[i]
            i += 1
    with open(dst, "w") as f:
        f.write(res)

