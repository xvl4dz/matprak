SetFactory("OpenCASCADE");

R = 1.5; // Major radius
r = 0.25; // Minor radius
lc = 0.1; // Mesh size

Torus(1) = {0,0,0, R, r, 2*Pi};

Surface Loop(1) = {1};
Volume(1) = {1};

// Mesh.Algorithm = 6; // Frontal
// Mesh.Algorithm3D = 4; // Delaunay
Mesh.CharacteristicLengthMax = lc;
Mesh.CharacteristicLengthMin = lc;
