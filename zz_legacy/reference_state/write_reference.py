import numpy as np

def write_reference(filename, radius, rho, dlnrho, d2lnrho, pressure,\
        temperature, dlnT, dsdr, entropy, gravity):
    f = open(filename, "wb")
    # may need to specify the data type for a successful read on
    # Rayleigh's end
    sigpi = np.array(314, dtype=np.int32)
    nr = np.array(len(radius), dtype=np.int32)
    f.write(sigpi.tobytes())
    # may need to specify the data type for a successful read on
    # Rayleigh's end
    sigpi = np.array(314, dtype=np.int32)
    nr = np.array(nr, dtype=np.int32)
    f.write(nr.tobytes())
    f.write(radius.tobytes())
    f.write(rho.tobytes())
    f.write(dlnrho.tobytes())
    f.write(d2lnrho.tobytes())
    f.write(pressure.tobytes())
    f.write(temperature.tobytes())
    f.write(dlnT.tobytes())
    f.write(dsdr.tobytes())
    f.write(entropy.tobytes())
    f.write(gravity.tobytes())
    f.close()
