# Package

version       = "0.1.0"
author        = "Benjamin Lee"
description   = "Homology-free viroid identification"
license       = "MIT"
srcDir        = "src"
installExt    = @["nim"]
bin           = @["viroid_search"]



# Dependencies

requires "nim >= 1.2.0"
requires "argpase >= 0.10.1"

task test, "Runs the test suite":
  exec "nim c -r -d:danger tests/test1"