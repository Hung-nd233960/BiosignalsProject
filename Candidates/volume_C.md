# Volume C: High-Resolution Time-Frequency Analysis of Real EEG Data

## C.1 First Look and Triage *(decide the direction from the data)*

### C.1.1 The EEG Dataset and Setup
To validate the signal-processing methodologies derived in Volume A and simulated in Volume B, we apply them to a real-world biomedical recording: a multi-channel infant electroencephalogram (EEG) stored in the standard European Data Format (`EEG.edf`). 

The recording contains 22 channels in total. The primary sampling frequency is $f_s = 200$ Hz. Initial channel inspection reveals that channels 25+, 26+, and 27+ are non-brain electrodes (analog auxiliary channels or external measurement devices) that do not conform to the International 10-20 System. The standard 10-20 scalp electrodes included in the recording are:
`Fp1`, `Fp2`, `F3`, `F4`, `C3`, `C4`, `P3`, `P4`, `O1`, `O2`, `F7`, `F8`, `T3`, `T4`, `T5`, `T6`, `Fz`, `Cz`, `Pz`.

### C.1.2 Time-Domain Visualization
The first stage of our pipeline involves loading the raw recording, stripping the `"EEG "` prefix from the channel names for readability, and plotting the first 100 seconds of the data. 

At this stage, the raw signals display substantial baseline drift (DC offset), low-frequency breathing artifacts, and high-frequency muscle (EMG) noise. Some channels exhibit brief, high-amplitude voltage shifts corresponding to physical movements or electrode contact adjustments.

### C.1.3 Frequency Triage and Dominant Bands
A preliminary Fast Fourier Transform (FFT) reveals that the spectral energy is distributed across the standard EEG bands:
- **Delta ($\delta$):** 1.0 – 4.0 Hz
- **Theta ($\theta$):** 4.0 – 8.0 Hz
- **Alpha ($\alpha$):** 8.0 – 13.0 Hz
- **Beta ($\beta$):** 13.0 – 30.0 Hz

For infants, the resting and sleep EEG is heavily dominated by low-frequency activity, particularly in the delta and theta bands, with intermittent emerging alpha rhythms (8-10 Hz) in the occipital region during quiet states.

### C.1.4 Triage Decision: Archetype Mapping
Based on this initial inspection, the infant EEG signal does not resemble a simple stationary single tone. Instead, it is characterized as:
- A mixture of slow background oscillations (delta/theta bands) and transient bursting rhythms (alpha band), mapped to the **Multi-Chirp (AA.4)** and **Mixed Tones (AA.2)** archetypes.
- Clicks, eye-blinks, and electrode pops, mapped to the **Transient/Pulse (AA.5)** archetype.
- Biological background activity and amplifier noise, mapped to the **Noise (AA.6)** archetype.

This combination of simultaneous frequencies, rapid temporal transitions, and noise makes the EEG dataset an ideal, highly challenging testbed for the high-resolution SPWVD method.

---

## C.2 Stationary Characterization — DFT and Band Power *(tone / mixed-tone)*

### C.2.1 Quantitative Pipeline
To establish a stationary baseline, we implement a quantitative band-power analysis using Welch's method (Section A.4.4). The pipeline consists of:
1. Creating a copy of the dataset for band-power quantification.
2. Applying a **Notch Filter at 50.0 Hz** to suppress power-line interference.
3. Applying a **Bandpass Filter from 1.0 Hz to 40.0 Hz** (FIR filter, Hamming window) to isolate the physiological range and eliminate baseline drift (DC) and high-frequency EMG noise.
4. Computing the Power Spectral Density (PSD) using Welch's method with a segment length of 5.0 seconds (frequency resolution $\Delta f = 0.2$ Hz) and a 50% overlap.
5. Integrating the PSD across the four frequency bands using the trapezoidal rule (`np.trapezoid`) to compute the absolute power ($\mu\text{V}^2$) and relative power (%) for each electrode.

### C.2.2 Results: Spatial Energy Distribution
We visualize the relative power distribution across all standard electrodes using a stacked bar chart. 

The analysis reveals that **Delta power is dominant across all channels**, accounting for 50% to 75% of the total spectral energy. This delta dominance is a typical neurophysiological signature of the infant brain. The occipital channels (`O1` and `O2`) show a clear, localized elevation in **Alpha relative power** (15% to 20%), indicating the presence of an active visual occipital alpha rhythm. 

---

## C.3 Time-Varying Characterization — STFT Spectrogram *(chirp / non-stationary)*

### C.3.1 Spectrogram Configuration
To observe how these rhythms evolve over time, we compute the Short-Time Fourier Transform (STFT) spectrogram of channel `O1` (the occipital representative). 

We select a window size of $\tau = 4.0$ s (800 samples at $f_s = 200$ Hz) with a 50% overlap. This window size is chosen specifically to balance frequency resolution (necessary to separate delta from theta bands) and time resolution (necessary to track the onset and termination of alpha bursts).

### C.3.2 Window Performance: Hanning vs. Hamming
We compare the spectrogram computed using a Hanning window against one using a Hamming window:
- Both window functions successfully suppress spectral leakage (Section A.3.4), replacing the sinc-like leakage ripples with smooth, localized energy bands.
- The **Hanning window** provides slightly superior side-lobe suppression, yielding a cleaner background spectrogram and making weak beta activity visible.
- The **Hamming window** offers a slightly narrower main lobe, resulting in marginally thinner horizontal frequency bands.

The spectrogram reveals that the occipital alpha rhythm is not continuous; instead, it appears as discrete **transient bursts** lasting between 3 and 10 seconds, coexisting with a continuous, low-frequency delta background.

---

## C.4 Events and Artifacts — Statistics and Transient Detection *(transient / noise)*

### C.4.1 Preprocessing Decisions
EEG recordings are highly sensitive to biological and environmental contamination:
1. **Drop-channel Triage:** We explicitly discard the auxiliary channels `25+`, `26+`, and `27+` at the start of the pipeline. Because these channels record non-brain potentials (ECG, breathing, electrode impedance tests), keeping them would distort the average reference.
2. **DC Baseline Filtering:** For Fast Fourier Transforms, baseline drift acts as a massive DC component. If left unmitigated, the leakage from the 0 Hz bin would spill into the delta band, inflating delta power estimates. We apply a sharp Highpass filter at 1.0 Hz to eliminate baseline drift, which keeps the 0 Hz energy spike at zero and ensures the integrity of physiological band power.

### C.4.2 Detection of Transient Spikes
Transient artifacts (such as eye blinks or motion clicks) appear in the raw time series as sudden, high-amplitude spikes. In the spectrogram, these map to vertical broadband stripes. By applying the statistical detection thresholds derived in Section A.4.2, we flag samples that exceed $3\sigma$ of the background noise floor, allowing these segments to be cleaned or excluded from stationary Welch averaging.

---

## C.5 High-Resolution Time-Frequency — WVD / SPWVD *(the payoff)*

### C.5.1 Segment Selection
The raw Wigner-Ville Distribution (WVD) cannot be applied to the entire, unedited multi-channel EEG recording. Because the full recording contains multiple simultaneous bands, muscle artifacts, and noise, the WVD's quadratic structure would generate a dense, uninterpretable "cross-term soup" (Section A.7.3).

To perform a meaningful high-resolution analysis, we select a **cleaned 2.0-second segment of channel `O1`** (from $t = 100.0$ s to $t = 102.0$ s) containing active occipital alpha waves. The signal is preprocessed using notch filtering (50 Hz), bandpass filtering (1–40 Hz), and average referencing. We then convert it to its analytic form and interpolate by 2 to prevent aliasing.

We compute and compare three distributions:
1. **STFT Spectrogram:** Hanning window of length 128 samples (0.64 s), overlapping by 50%.
2. **Raw WVD:** Calculated using the analytic signal with no time or frequency smoothing.
3. **SPWVD:** Calculated using a frequency lag window $h$ (Hann, length 51) and a time smoothing window $g$ (Hann, length 15).

### C.5.2 Results and Visual Analysis
We implement a **Dual-Stack Visualization** of the results, plotting the time-frequency distributions on a linear scale in Figure C.1 and on a logarithmic decibel (dB) scale in Figure C.2.

![Figure C.1 - Real EEG Time-Frequency Analysis (Linear Scale)](../results/graphs/volume_c/figure_C_01.png)

![Figure C.2 - Real EEG Time-Frequency Analysis (dB Scale)](../results/graphs/volume_c/figure_C_02.png)

**1. STFT Spectrogram (Top Panel):**
- Displays the alpha rhythm as a single, thick, blurred band centered around 10–11 Hz.
- The temporal changes in the alpha amplitude are smeared across the columns, and the exact frequency shifts are difficult to resolve.
- The low-frequency delta background is merged with the alpha band in a low-resolution cloud.

**2. Raw WVD (Middle Panel):**
- The trajectories of the components are extremely sharp, but the time-frequency plane is completely unreadable.
- A dense web of time-oscillating and frequency-oscillating cross-term ghosts fills the space between 0 Hz and 20 Hz, completely masking the actual biological rhythms.
- On the dB scale (Figure C.2, middle panel), these cross-terms cover the entire dynamic range, illustrating the collapse of the raw WVD in multi-component scenarios.

**3. SPWVD (Bottom Panel):**
- By applying the lag window $h=51$ and the time window $g=15$, the SPWVD successfully **suppresses the cross-term ghosts**.
- It reveals **two distinct, sharp trajectories**: a slow, background delta rhythm under 4 Hz, and an active occipital alpha rhythm oscillating between 10 Hz and 12 Hz.
- Unlike the blurred STFT, the SPWVD displays the alpha rhythm as a highly localized, undulating wave. The frequency modulations (drifting between 10.5 Hz and 11.5 Hz) and the precise onset/offset timing of the rhythm are visible.
- The dB scale (Figure C.2, bottom panel) confirms that the noise floor and ghost artifacts have been suppressed below $-30$ dB, while the high concentration of the true paths is preserved.

### C.5.3 Tradeoff Characterization
This analysis illustrates that the SPWVD successfully delivers on its theoretical promise. By decoupling the time and frequency smoothing parameters, it allows us to suppress cross-term artifacts while maintaining a level of time-frequency concentration that is impossible to achieve with the STFT. 

However, this resolution comes at the cost of computational complexity and requires careful, manual tuning of the windows $g$ and $h$ based on the frequency separation of the expected rhythms.

---

## C.6 Synthesis — What We Found

### C.6.1 Methodological Match
Our work on real infant EEG data confirms that different DSP tools are suited to different neurophysiological structures, validating the taxonomy developed in Appendix AA:
- **Welch's PSD (C.2):** Best suited for long-term, global stationary characterization (e.g., computing relative band power tables across the scalp).
- **STFT Spectrogram (C.3):** Highly robust and computationally efficient for broad time-varying characterization, but limited by Heisenberg blurring.
- **SPWVD (C.5):** The optimal choice for high-resolution analysis of selected, short segments containing overlapping or frequency-modulating rhythms (such as occipital alpha bursts or sleep spindles). It tracks trajectories with exceptional clarity once cross-terms are managed.

### C.6.2 Defensible Claims
We conclude that:
1. The infant EEG dataset studied is delta-dominated with a clear occipital alpha rhythm localized in space (`O1`/`O2` channels) and time (bursting behavior).
2. Naive WVD is unusable on real EEG due to cross-term interference, but the Hilbert analytic signal and SPWVD smoothing provide a highly effective, high-resolution solution.
3. The choice of time-frequency transform must be driven by the specific clinical or engineering question: global power distribution calls for Welch; trajectory tracking of short rhythms calls for SPWVD.

---

## Appendix C — EEG Background

### C.1 EEG Signals and Montages
An Electroencephalogram (EEG) measures the electrical activity of the brain via electrodes placed on the scalp. The microvolt-scale potentials originate from the postsynaptic currents of thousands of aligned cortical pyramidal neurons. 

The standard positioning of these electrodes is governed by the **International 10-20 System**, which uses proportional anatomical distances (10% and 20%) between the nasion, inion, and preauricular points to ensure consistent placement:
- **Fp (Frontopolar):** Forehead region.
- **F (Frontal):** Frontal lobe, associated with motor planning and cognitive function.
- **C (Central):** Motor cortex and somatosensory strip.
- **P (Parietal):** Parietal lobe, associated with sensory integration.
- **O (Occipital):** Occipital lobe, the primary visual processing center.
- **T (Temporal):** Temporal lobe, associated with auditory processing and memory.

Odd-numbered electrodes (e.g., `Fp1`, `O1`) are placed on the left hemisphere, even-numbered electrodes (e.g., `Fp2`, `O2`) on the right hemisphere, and the letter `z` (e.g., `Fz`, `Cz`, `Pz`) denotes midline electrodes.

### C.2 Physiological Rhythms and States
Different frequency bands are associated with distinct cognitive and physiological states in the adult baseline:
- **Delta (0.5 – 4 Hz):** Prominent during deep, slow-wave sleep (Non-REM Stage 3). High amplitude delta in awake adults is typically pathological.
- **Theta (4 – 8 Hz):** Associated with drowsiness, light sleep, and memory consolidation.
- **Alpha (8 – 13 Hz):** The dominant rhythm of awake, relaxed adults. It is strongest in the occipital region when the eyes are closed and is immediately blocked (attenuated) when the eyes open or during cognitive tasks.
- **Beta (13 – 30 Hz):** Associated with active thinking, focus, and motor execution.
- **Gamma (>30 Hz):** Associated with high-level cognitive processing and sensory integration.

### C.3 Infant and Neonatal EEG Characteristics
The EEG of infants and neonates differs fundamentally from the adult standard:
1. **Delta Dominance:** The infant EEG is dominated by high-amplitude slow-wave activity (delta and theta), reflecting the ongoing development and myelination of cortical pathways.
2. **Discontinuity:** Sleep EEG in premature infants and neonates exhibits alternating periods of high-voltage activity and low-voltage quiescent periods, a normal physiological pattern known as **discontinuous EEG** (e.g., *tracé alternant* or *tracé discontinu*).
3. **Maturity Marker:** The presence and characteristics of these discontinuous patterns serve as a key marker of brain maturity and gestational age. As the infant matures, the quiet periods shorten, and the EEG transitions to a continuous pattern.

### C.4 Common Artifacts in EEG
EEG recordings are highly sensitive to biological and environmental contamination:
- **EOG (Electrooculogram):** Eye movements and blinks create massive, low-frequency voltage deflections, concentrated in the frontal channels (`Fp1`/`Fp2`).
- **EMG (Electromyogram):** Muscle tension in the jaw or neck generates high-frequency, broadband noise (>30 Hz) that can mask beta and gamma bands.
- **ECG (Electrocardiogram):** Heartbeat potentials can leak into the recording, appearing as sharp, periodic spikes across multiple electrodes.
- **Electrode Pop:** A sudden shift in contact impedance creates a sharp, high-amplitude transient that resembles a biological spike.

### C.5 Dataset Specifications
The dataset analyzed in this volume is a real EEG recording (`EEG.edf`) of an infant subject. The recording has a sampling rate of 200 Hz. The scalp electrodes are referenced to a computed average reference over the full electrode set, Notch-filtered at 50 Hz, and Bandpass-filtered at 1–40 Hz. The analysis in Section C.5 focuses on the `O1` channel, where occipital alpha bursts are most prominent.

### References
1. Niedermeyer, E., & da Silva, F. L. (2005). *Electroencephalography: Basic Principles, Clinical Applications, and Related Fields*. Lippincott Williams & Wilkins.
2. de Haan, M. (2007). *Infant EEG and Event-Related Potentials*. Psychology Press.
3. Scher, M. S. (2006). *Ontogeny of EEG Sleep Patterns in the Human Infant*. *Journal of Clinical Neurophysiology*, 23(3), 229-241.
