#!/usr/bin/env bash
#
# Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

#path is never empty
PATH="/home/yoar/votca/xtp/votca/bin:${PATH}"
#debian wants to have : at the end
MANPATH="/home/yoar/votca/xtp/votca/share/man:${MANPATH}"
LD_LIBRARY_PATH="/home/yoar/votca/xtp/votca/lib${LD_LIBRARY_PATH:+:}${LD_LIBRARY_PATH}"

VOTCASHARE="/home/yoar/votca/xtp/votca/share/votca"

export PATH MANPATH LD_LIBRARY_PATH VOTCASHARE

#votca rc files
for rc in "${VOTCASHARE}"/rc/*rc.bash; do
  [ -r "$rc" ] && . "$rc"
done
unset rc

#bash completion
if [ -n "$BASH_VERSION" ]; then 
  for comp in "${VOTCASHARE}"/rc/*completion.bash; do
    [ -r "$comp" ] && source "$comp"
  done
  unset comp
fi
