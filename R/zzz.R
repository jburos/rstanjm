# Part of the rstanjm package
# Copyright (C) 2015, 2016 Trustees of Columbia University
# Copyright (C) 2016 Sam Brilleman
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

.onLoad <- function(libname, pkgname) { # nocov start
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) Rcpp::loadModule(m, what = TRUE)
} # nocov end

.onAttach <- function(...) {
  rstanjmLib <- dirname(system.file(package = "rstanjm"))
  pkgdesc <- utils::packageDescription("rstanjm", lib.loc = rstanjmLib)
  builddate <- gsub(';.*$', '', pkgdesc$Packaged)
  #packageStartupMessage(paste0("rstanjm (Version ", pkgdesc$Version, ")"))
  packageStartupMessage("* By default rstanjm fits models using a single MCMC chain. If you have a")
  packageStartupMessage("  multicore CPU with excess RAM, you can fit multiple MCMC chains in parallel")
  packageStartupMessage("  by specifying the 'chains' and 'cores' arguments when fitting your model.")
}

