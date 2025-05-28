class LightConeParticle {
    public:
        double x, y, z, vx, vy, vz;
};
extern std::vector<LightConeParticle> shell_particles;

void writeBinaryVTK(const std::string& filename,
                    const std::vector<LightConeParticle>& shell_particles);
void writeBinarySimple(const std::string& filename,
                    const std::vector<LightConeParticle>& shell_particles);
