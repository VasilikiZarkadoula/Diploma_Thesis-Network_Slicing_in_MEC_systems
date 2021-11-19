function h = randomGain(channels)
rng shuffle;

% exponential channel gain corresponding to Rayleigh fading
h = exprnd(2,channels,channels);
  