/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_CHECKPOINT_H
#define VOTCA_XTP_CHECKPOINT_H

#include <H5Cpp.h>
#include <votca/xtp/checkpoint_utils.h>
#include <votca/xtp/checkpointreader.h>
#include <votca/xtp/checkpointtable.h>
#include <votca/xtp/checkpointwriter.h>

namespace votca {
namespace xtp {

enum class CheckpointAccessLevel {
  READ = 0,    // only read no write access
  MODIFY = 1,  // if file exists, change it
  CREATE = 2   // create new file
};

std::ostream& operator<<(std::ostream& s, CheckpointAccessLevel l);

class CheckpointFile {
 public:
  CheckpointFile(std::string fileName);
  CheckpointFile(std::string fileName, CheckpointAccessLevel access);

  std::string getFileName();
  std::string getVersion();

  H5::H5File getHandle();

  CheckpointWriter getWriter();
  CheckpointWriter getWriter(const std::string _path);
  CheckpointReader getReader();
  CheckpointReader getReader(const std::string _path);

 private:
  std::string _fileName;
  H5::H5File _fileHandle;
  CptLoc _rootLoc;
  CheckpointAccessLevel _accessLevel;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_CHECKPOINT_H
