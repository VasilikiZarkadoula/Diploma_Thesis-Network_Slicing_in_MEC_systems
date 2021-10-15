function h_mMTC = randomGain(channels)

rng shuffle;
h_mMTC = exprnd(2,channels,channels);
  