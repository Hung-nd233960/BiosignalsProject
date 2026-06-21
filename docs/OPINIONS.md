# Some opinions that requires pondering and verdict

- [x] Lab B.3 should specify more on how scipy can find the maxima.
  - Resolved: added comment explaining `argrelextrema` (compares each point to `order` neighbors on each side) in volume_B.md.

- [x] In all reports, saying common modes like "fit a line" is a bit shallow. Specify which Linear Regression mode we use. and what function is it, if it is external, maybe brief on how they work.
  - Resolved: added OLS explanation for `linregress` in volume_B.md (Lab 3 decay regression) and volume_C.md (C.2.2 1/f fit). CLAUDE.md rule updated: "When using external functions, state briefly what they do and what method they use."
