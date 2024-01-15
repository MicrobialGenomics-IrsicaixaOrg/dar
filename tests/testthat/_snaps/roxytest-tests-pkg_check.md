# Function recipes_pkg_check() @ L18

    Code
      dar:::recipes_pkg_check(dar:::required_pkgs_error(), "step_aldex()")
    Message
      i 2 packages are needed for step_aldex() and are not installed: (randompackage, packrandom)
      * Start a clean R session then run: BiocManager::install(c("randompackage", "packrandom"))

---

    Code
      dar:::recipes_pkg_check(dar:::required_pkgs_aldex(), "step_aldex()")

