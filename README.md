# Eris creativity model

This project uses [Eris](https://github.com/erisproject/eris) to model the
creative process, including creation, production, and purchases of copies of a
creative work.

This is primarily modelled as authors, publishers, and readers, but that
terminology is for convenience: the model should be applicable to any sort of
creative output where copies of an original are the marketable good.

## Requirements

- [Eris](https://github.com/erisproject/eris)
- [gtkmm-3.0](http://www.gtkmm.org) (optional--required for GUI interface)
- A C++ compiler supporting the C++11 standard, such as
  [clang](http://clang.llvm.org/) (3.3+) or [g++](https://gcc.gnu.org/) (4.9+)
- [Boost](http://www.boost.org), including the filesystem library component.
- [TCLAP](http://tclap.sf.net) command-line parsing library.

## Compiling

To compile on a unix-like system, create a new build directory somewhere, then
from this directory run:

    cmake /path/to/creativity
    make -j4

To compile without the GTK GUI, change the cmake line to:

    cmake -DGUI=off /path/to/creativity

You can install directly to the system (usually under /usr/local) using:

    make install

or alternatively simply run the executables directly from the build directory.
To build a .deb to install on a Debian derivative, run:

    make package

followed by an appropriate package command to install the package, such as:

    dpkg -i eris_x.y.z~gityyyymmdd~tag_amd64.deb

You may also be able to use:

    cpack -G RPM

to generate a .rpm package instead of a .deb, but this functionality is
untested by the author.

## License

Copyright © 2015 Jason Rhinelander

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
