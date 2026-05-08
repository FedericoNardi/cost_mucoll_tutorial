# Muon Collider Data Generation Tutorial

*F. Nardi — COST Tutorial, May 2026*

---

## The Muon Collider concept

The Muon Collider (MC) is a proposed circular lepton collider with several key advantages:

- Less synchrotron radiation compared to electron colliders
- No hadronic pile-up
- Able to reach multi-TeV centre-of-mass energies (current performance studies target √s = 10 TeV)

Proposed sites are CERN and Fermilab, with a collider ring of O(10 km).

### Challenges

**Beamline cooling** — the cooling stage must produce a coherent muon beam to feed the acceleration chain. 6D-magnet cooling is the most promising candidate.

**Beam-Induced Background (BIB)** — muons are unstable (τ ≈ 2.2 μs, boosted to ~2.2 ms at 10 TeV) and decay along the beamline, producing a cloud of secondary particles (mostly photons and neutrons) that interact with the detector. A double-cone tungsten nozzle at the Machine-Detector Interface screens most of the radiation, but the interaction point remains unshielded.

---

## Software stack

This tutorial follows the [MuCol software wiki](https://mcd-wiki.web.cern.ch/software/tutorials/) and is based on the **ILCSoft** and **Key4hep** software suites. It requires a container manager (e.g. [Apptainer](https://apptainer.org/docs/user/latest/)) and access to CVMFS.

### Get the steering files

```bash
git clone https://github.com/FedericoNardi/cost_mucoll_tutorial.git
```

### Start the container and source the environment

```bash
apptainer run -B <disk/mount/location/> \
  /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/muoncollider/mucoll-deploy/mucoll:2.9-alma9

source /opt/setup_mucoll.sh
```

The `-B` flag mounts additional storage from your host machine into the container.

---

## Part 1 — Particle Gun

A particle gun generates single-particle events with user-defined kinematics. This is the simplest way to test the simulation chain.

### 1.1 Event generation

Generate a `.slcio` file containing initialised MCParticles:

```bash
python mucoll-benchmarks/generation/pgun/pgun_lcio.py \
  --seed 42 \        # random seed
  -e 10 \            # number of events
  --pdg 11 \         # particle type (11 = electron)
  --p 100 \          # momentum [GeV]
  --theta 10 170 \   # polar angle range [deg]
  -- output_gen.slcio
```

Inspect the output:

```bash
anajob output_gen.slcio
```

### 1.2 Detector simulation

Run Geant4 simulation via `ddsim`. The `MuColl_v1` geometry for the 3TeV detector studies is available at `$MUCOLL_GEO`:

```bash
echo $MUCOLL_GEO
# /opt/.../lcgeo/.../MuColl_v1/MuColl_v1.xml
```

```bash
ddsim --steeringFile mucoll-benchmarks/simulation/ilcsoft/steer_baseline.py \
  --inputFile  output_gen.slcio \
  --outputFile output_sim.slcio
```

The output file contains:

- `MCParticle` — input and secondary particles produced during detector interaction
- `SimTrackerHit` / `SimCalorimeterHit` — energy deposits per particle

```bash
anajob output_sim.slcio
```

### 1.3 Digitisation

Convert SimHits into realistic detector hits, applying energy thresholds, smearing, and timing cuts. Configuration lives in `digi_steer.py` (e.g. ECal threshold `5e-5 GeV`, timing window `[-0.5, 15] ns`).

```bash
k4run mucoll-benchmarks/digitisation/k4run/digi_steer.py \
  --LcioEvent.Files <path/for/sim/file>/output_sim.slcio \
  --outputPath <path/for/digitised/file/>
```

Output: `output_digi.slcio` (full collections) and `output_digi_light.slcio` (reduced).

Downstream reconstruction (tracks, clusters, PFOs) uses ACTS for tracking and Pandora for particle flow. See the [MuCol wiki](https://mcd-wiki.web.cern.ch) for details.

---

## Part 2 — Monte Carlo signal events

Instead of particle guns, realistic physics events can be generated with a Monte Carlo generator and passed to `ddsim` as HepMC3 or `stdhep` input.

### 2.1 Signal process

The benchmark signal is:

```
μ⁺μ⁻ → H νν → bb̄ νν    (WW-fusion, √s = 10 TeV)
```

WHIZARD handles matrix element integration and phase space sampling; Pythia8 performs parton showering and hadronisation.

### 2.2 Event generation with WHIZARD

Run in a dedicated folder to keep the compilation files isolated:

```bash
mkdir H_bb && cd H_bb
whizard ../mucoll-benchmarks/generation/signal/whizard/mumu_H_bb_3TeV.sin
```

Output: `mumu_H_bb_3TeV.hepmc`

**Integration settings** — the `.sin` file controls the phase space grid optimisation:

```
integrate (hdec) { iterations = 3:1000:"gw", 1:5000 }   # tutorial (fast, ~seconds)
integrate (hdec) { iterations = 20:8000000:"gw", 5:20000000 }  # production (hours)
```

For a tutorial with ~10 events, the light settings are sufficient. 
Alternatively, you can download a cached .grid file from 
```bash
cp /cvmfs/muoncollider.cern.ch/datasets/tutorial_20230705/gen_Hbb/mumu_H_bb_3TeV.hepmc ./
```
or (if `/cvmfs` is missing)
```bash
wget https://nbartosi.web.cern.ch/tutorial_20230705/gen_Hbb/mumu_H_bb_3TeV.hepmc
```
and generate the MC particles without compiling matrix elements:
```bash
whizard --no-integration mumu_H_bb_3TeV.sin
````

### 2.3 Detector simulation and digitisation

Same commands as the particle gun case, pointing to the HepMC3 file:

```bash
ddsim --steeringFile mucoll-benchmarks/simulation/ilcsoft/steer_baseline.py \
  --inputFile  H_bb/mumu_H_bb_3TeV.hepmc \
  --outputFile output_sim.slcio

k4run mucoll-benchmarks/digitisation/k4run/digi_steer.py \
  --LcioEvent.Files output_sim.slcio
```

---

## Part 3 — BIB overlay

BIB particles are simulated separately (FLUKA → `ddsim`) and stored as pre-simulated SimHit files. They are overlaid onto the signal SimHits during digitisation.

### Overlay processor

`OverlayTimingRandomMix` draws independently from μ⁺ and μ⁻ BIB file pools and randomly composes each BIB event on the fly, providing genuine event-to-event stochasticity even with a small file pool. Both arguments expect a **directory** of `.slcio` files.

```bash
export MUPLUS="/path/to/bib/sim_mu_plus/"
export MUMINUS="/path/to/bib/sim_mu_minus/"

k4run mucoll-benchmarks/digitisation/k4run/digi_steer.py \
  --LcioEvent.Files            output_sim.slcio \
  --doOverlayFull \
  --OverlayFullPathToMuPlus    $MUPLUS \
  --OverlayFullPathToMuMinus   $MUMINUS \
  --OverlayFullNumberBackground 192     # 192 = full bunch crossing (45 φ-clones/particle)
```

### Timing cuts

BIB hits are filtered per subdetector by arrival time. Only hits within the integration window are retained — this is the primary BIB rejection mechanism at hit level.

| Subdetector | Window |
|---|---|
| Vertex / inner tracker | ±0.2–0.5 ns |
| Outer tracker | ±0.5 ns |
| Calorimeters | −0.5 to +15 ns |

---

## File format reference

| Extension | Format | Tool |
|---|---|---|
| `.slcio` | LCIO (legacy) | `anajob`, `lcio_event_counter` |
| `.hepmc` | HepMC3 | WHIZARD output, `ddsim` input |
| `.edm4hep.root` | EDM4hep (Key4hep native) | `k4run` output |

---

## References

- [MuCol software wiki — CERN 2023 tutorial](https://mcd-wiki.web.cern.ch/software/tutorials/cern2023)
- [MuCol software wiki — Fermilab 2024 tutorial](https://mcd-wiki.web.cern.ch/software/tutorials/fermilab2024)
- [mucoll-benchmarks repository](https://github.com/MuonColliderSoft/mucoll-benchmarks)
- [Key4hep](https://key4hep.github.io/key4hep-doc/)
- [ILCSoft](https://github.com/iLCSoft)