#pragma once
#include "AtomData.h"
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace DistTool {

/**
 * Reads LAMMPS "dump" files (text format).
 * Supports: id type x y z  (also xu yu zu for unwrapped coords).
 * Can iterate frame-by-frame or read all frames at once.
 */
class LammpsDumpReader {
public:
    explicit LammpsDumpReader(const std::string& path);

    // Returns true while unread frames remain.
    bool hasNext() const;

    // Read and return the next frame.  Throws std::runtime_error at EOF.
    Frame readNext();

    // Convenience: read every frame in the file.
    std::vector<Frame> readAll();

    // Thrown by parseFrame() only when EOF is reached before finding any
    // ITEM: TIMESTEP header — i.e. a clean end-of-file between frames.
    // readAll() catches this to terminate iteration; readNext() lets it propagate
    // as a std::runtime_error so callers see a meaningful message.
    struct CleanEof : std::runtime_error {
        CleanEof() : std::runtime_error("No more frames in file") {}
    };

private:
    std::string           path_;
    mutable std::ifstream file_;  // mutable: peek() in hasNext() updates stream state

    Frame parseFrame();

    // Parse the column list after "ITEM: ATOMS"
    static std::vector<std::string> splitColumns(const std::string& header);
};

} // namespace DistTool
