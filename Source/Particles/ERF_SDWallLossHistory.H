#ifndef ERF_SD_WALL_LOSS_HISTORY_H_
#define ERF_SD_WALL_LOSS_HISTORY_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <memory>
#include <string>

namespace sdwl {

/** Host-only history for one wall-loss-enabled particle container. All ranks
 * call in the same order, but only the IOProcessor owns state or accesses files.
 * Nothing in this class feeds back into the particles; no broadcast is needed.
 * Pass a stable, single-line signature identifying the physics and bin edges.
 * Writers serialize validation and output with a short-lived .lock directory.
 * After an interrupted write, check the history and that no writer is active
 * before removing a leftover lock. A lock is never stolen automatically.
 */
class History
{
public:
    History ();
    ~History ();
    History (const History&) = delete;
    History& operator= (const History&) = delete;

    /** Once per accepted step, starting at 1. time is the START of that step.
     * global on the writer is already reduced: side/top/bottom number, then
     * total material mass [kg] in the same order, then 1..64 bin counts.
     * Every step contributes mass. Positive multiples of cadence emit a row;
     * cadence zero suppresses text only. Counts/rates are not interval means.
     */
    void record (const std::string& name, const std::string& signature, int cadence,
                 int step, amrex::Real time, amrex::Real dt,
                 const amrex::Vector<amrex::Real>& global);

    /** After particle Checkpoint has created dir/name. Legal before record(1):
     * the lineage is established now; cadence and bin count are not yet known.
     * Publishes an immutable sidecar; writing identical bytes again is a no-op.
     * Does not consult the native ERF Header during checkpoint construction.
     */
    void checkpoint (const std::string& dir, const std::string& name,
                     const std::string& signature) const;

    /** Requires a sidecar even at cadence zero. With text enabled, verifies the
     * exact expected segment or exclusively creates a marked continuation in
     * the working directory. Never replays the checkpoint's sampled row.
     * If present, the native Header supplies level-zero accepted step/end time;
     * without it a step-zero particle-only checkpoint has no initial clock.
     */
    void restart (const std::string& dir, const std::string& name,
                  const std::string& signature, int cadence);

private:
    struct State;
    State& bind (const std::string& name, const std::string& signature) const;
    // Lazy writer identity: even two pre-step checkpoints must share a lineage.
    mutable std::unique_ptr<State> m_state;
};

} // namespace sdwl

#endif
