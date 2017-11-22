from ZZMatrixElement.MELA.mela import Mela, SimpleParticleCollection_t

class LHEFile_JHUGenVBFVH(object):
    """
    Simple class to iterate through an LHE file and calculate probabilities for each event
    For JHUGen VBF and VH only
    """
    mela = {}
    def __init__(self, filename, *melaargs):
        self.filename = filename
        if melaargs not in type(self).mela:
            type(self).mela[melaargs] = Mela(*melaargs)
        self.mela = type(self).mela[melaargs]
        self.f = open(self.filename)
    def __enter__(self, *args, **kwargs):
        self.f.__enter__(*args, **kwargs)
        return self
    def __exit__(self, *args, **kwargs):
        return self.f.__exit__(*args, **kwargs)

    def __iter__(self):
        event = ""
        for linenumber, line in enumerate(self.f, start=1):
            if "<event>" not in line and not event:
                continue
            event += line
            if "</event>" in line:
                try:
                    particles = event.strip().split("\n")[2:-1]  #remove the first 2 lines, <event> and event info, and the last line, </event>
                    self.daughters = SimpleParticleCollection_t([particle for particle in particles if int(particle.split()[0]) == 25])
                    self.associated = SimpleParticleCollection_t([particle for particle in particles
                                                                   if abs(int(particle.split()[0])) in (1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16)
                                                                   and int(particle.split()[1]) == 1])
                    self.mothers = SimpleParticleCollection_t([particle for particle in particles if int(particle.split()[1]) == -1])
                    assert len(self.daughters) == 1 and len(self.associated) == len(self.mothers) == 2
                    self.mela.setInputEvent(self.daughters, self.associated, self.mothers, True)
                    yield self
                    event = ""
                except GeneratorExit:
                    raise
                except:
                    print "On line", linenumber
                    raise
                finally:
                    try:
                        self.mela.resetInputEvent()
                    except:
                        pass

    def __getattr__(self, attr):
        return getattr(self.mela, attr)
    def __setattr__(self, attr, value):
        if attr in ("filename", "f", "mela", "daughters", "mothers", "associated"):
            super(LHEFile_JHUGenVBFVH, self).__setattr__(attr, value)
        else:
            setattr(self.mela, attr, value)
