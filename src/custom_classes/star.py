from amuse.lab import units, constants


class Star:
    def __init__(self, 
                particles, 
                timestep, 
                stellar_evo_code=None):
        self.particles = particles
        self.model_time = 0 | units.yr
        self.timestep = timestep

        self.stellar = stellar_evo_code
        if self.stellar:
            self.stellar.particles.add_particles(self.particles)
            self.channel_mass = self.stellar.particles.new_channel_to(self.particles)

    def pre_evolve(self, end_time):
        if not self.stellar:
            print("This star does not have a stellar evo code!")
            return
        self.stellar.evolve_model(end_time)
        self.channel_mass.copy()

    def get_gravity_at_point(self, eps, x, y, z): 
        x -= self.particles.x
        y -= self.particles.y
        z -= self.particles.z
        
        r = (x**2+y**2+z**2).sqrt()
        a = - constants.G * self.particles.mass/r**2
          
        ax = a*x/r
        ay = a*y/r
        az = a*z/r
        
        return ax, ay, az

    def update_mass(self):
        if self.stellar:
            self.stellar.evolve_model(self.stellar.model_time + self.timestep)
            self.channel_mass.copy()

    def update_pos(self):
        self.particles.x += self.particles.vx * self.timestep
        self.particles.y += self.particles.vy * self.timestep
        self.particles.z += self.particles.vz * self.timestep

    def evolve_model(self, time):
        while self.model_time < time:
            self.update_pos()
            self.update_mass()
            self.model_time += self.timestep
        

