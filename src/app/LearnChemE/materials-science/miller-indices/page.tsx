'use client';

import React, { useState, useMemo, Suspense } from 'react';
import * as THREE from 'three';
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Slider } from "@/components/ui/slider";
import { Label } from "@/components/ui/label";
import { Checkbox } from "@/components/ui/checkbox";
import { Input } from "@/components/ui/input";
import { Canvas } from '@react-three/fiber';
import { OrbitControls, Html } from '@react-three/drei';

const latticeConstant = 1.0; // Define a unit for lattice constant

interface Atom {
  x: number;
  y: number;
  z: number;
  radius: number;
  color: string;
}

interface CrystalSceneProps {
  crystalStructure: string;
  repeats: number;
  radiusPrimary: number;
  radiusSecondary: number;
  shiftX: number;
  shiftY: number;
  shiftZ: number;
  showUnitCell: boolean;
  showMillerPlane: boolean;
  millerH: string;
  millerK: string;
  millerL: string;
  unitCellSceneOrigin: THREE.Vector3; // Added prop
}

function CrystalScene({
  crystalStructure,
  repeats,
  radiusPrimary,
  // radiusSecondary, // Not actively used for BCC/FCC atom differentiation yet
  shiftX,
  shiftY,
  shiftZ,
  showUnitCell,
  showMillerPlane,
  millerH,
  millerK,
  millerL,
  unitCellSceneOrigin // Consumed as prop
}: CrystalSceneProps) {

  const atoms = useMemo(() => { // Changed: unitCellSceneOrigin removed from return
    const generatedAtoms: Atom[] = [];
    let basePositions: { x: number; y: number; z: number }[] = [];

    switch (crystalStructure) {
      case 'SC':
        basePositions = [{ x: 0, y: 0, z: 0 }];
        break;
      case 'BCC':
        basePositions = [
          { x: 0, y: 0, z: 0 },    // Corner
          { x: 0.5, y: 0.5, z: 0.5 } // Body center
        ];
        break;
      case 'FCC':
        basePositions = [
          { x: 0, y: 0, z: 0 },    // Corner
          { x: 0.5, y: 0.5, z: 0 }, // Face center XY
          { x: 0.5, y: 0, z: 0.5 }, // Face center XZ
          { x: 0, y: 0.5, z: 0.5 }  // Face center YZ
        ];
        break;
      default:
        basePositions = [{ x: 0, y: 0, z: 0 }];
    }

    // This offset centers the entire supercell block if repeats > 1.
    // If repeats = 1, offset is 0.
    const supercellCenteringOffset = (repeats > 1) ? (repeats * latticeConstant - latticeConstant) / 2 : 0;
    
    // unitCellSceneOrigin is now passed as a prop, so its direct calculation here is removed.
    // The atom positions will be relative to an ideal lattice origin (0,0,0)
    // and then the entire scene (including the unit cell box and miller plane)
    // will be effectively shifted by unitCellSceneOrigin in the parent component
    // or by how unitCellSceneOrigin is used in positioning those elements.
    // For atom generation, we use the supercellCenteringOffset to place them correctly
    // relative to the scene's center, assuming the unitCellSceneOrigin prop handles
    // the (0,0,0) unit cell's corner position.

    for (let rX = 0; rX < repeats; rX++) {
      for (let rY = 0; rY < repeats; rY++) {
        for (let rZ = 0; rZ < repeats; rZ++) {
          // Corner of the current cell in ideal lattice coordinates (before any centering)
          const cellCornerX_ideal = rX * latticeConstant;
          const cellCornerY_ideal = rY * latticeConstant;
          const cellCornerZ_ideal = rZ * latticeConstant;

          basePositions.forEach(bp => {
            generatedAtoms.push({
              // Atom position: ideal cell corner + (basis pos + fractional shift) * lc - supercell centering
              x: cellCornerX_ideal + (bp.x + shiftX) * latticeConstant - supercellCenteringOffset,
              y: cellCornerY_ideal + (bp.y + shiftY) * latticeConstant - supercellCenteringOffset,
              z: cellCornerZ_ideal + (bp.z + shiftZ) * latticeConstant - supercellCenteringOffset,
              radius: radiusPrimary, // For now, all atoms use primary radius
              color: '#FF6B6B' // Example: Red for primary
              // TODO: Could use radiusSecondary and different colors for different basis atoms if needed
            });
          });
        }
      }
    }
    return generatedAtoms; // Changed: only return atoms
  }, [crystalStructure, repeats, radiusPrimary, shiftX, shiftY, shiftZ]); // unitCellSceneOrigin removed from dependencies here

  const unitCellEdges = useMemo(() => {
    if (!showUnitCell) return null;
    
    const points: THREE.Vector3[] = [];
    const ucX = unitCellSceneOrigin.x; 
    const ucY = unitCellSceneOrigin.y; 
    const ucZ = unitCellSceneOrigin.z;

    const a = latticeConstant;
    const p = [ // Vertices of the unit cell
      new THREE.Vector3(ucX,     ucY,     ucZ), new THREE.Vector3(ucX + a, ucY,     ucZ), 
      new THREE.Vector3(ucX + a, ucY + a, ucZ), new THREE.Vector3(ucX,     ucY + a, ucZ),
      new THREE.Vector3(ucX,     ucY,     ucZ + a), new THREE.Vector3(ucX + a, ucY,     ucZ + a), 
      new THREE.Vector3(ucX + a, ucY + a, ucZ + a), new THREE.Vector3(ucX,     ucY + a, ucZ + a)
    ];

    // Edges
    [ p[0],p[1], p[1],p[2], p[2],p[3], p[3],p[0], // bottom face
      p[4],p[5], p[5],p[6], p[6],p[7], p[7],p[4], // top face
      p[0],p[4], p[1],p[5], p[2],p[6], p[3],p[7]  // vertical edges
    ].forEach(point => points.push(point));
    
    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    return geometry;

  }, [showUnitCell, unitCellSceneOrigin]); // unitCellSceneOrigin is now a prop

  const millerPlaneData = useMemo(() => {
    if (!showMillerPlane) return null;
    
    const h_val = parseInt(millerH);
    const k_val = parseInt(millerK);
    const l_val = parseInt(millerL);

    if (isNaN(h_val) || isNaN(k_val) || isNaN(l_val) || (h_val === 0 && k_val === 0 && l_val === 0)) {
      return null; 
    }

    const normalVec = new THREE.Vector3(h_val, k_val, l_val);
    if (normalVec.lengthSq() === 0) return null; // Should be caught by above, but defensive
    
    // Distance from the unit cell's origin to the plane (for N=1 family member)
    // d = a / |(h,k,l)| if plane eq is hx+ky+lz = a
    const distToPlaneFromCellOrigin = latticeConstant / normalVec.length(); 
    const planeNormalNormalized = normalVec.clone().normalize();
    
    // Position of the plane's center, relative to the visualized unit cell's origin
    const planeCenterRelativeToCellOrigin = planeNormalNormalized.clone().multiplyScalar(distToPlaneFromCellOrigin);

    // Final absolute position of the plane's center in the scene
    const finalPlaneCenterPosition = unitCellSceneOrigin.clone().add(planeCenterRelativeToCellOrigin);
    
    // Make plane large enough to be seen across the supercell
    const planeSize = repeats * latticeConstant * Math.max(1.5, Math.abs(h_val||1), Math.abs(k_val||1), Math.abs(l_val||1));
    const geometry = new THREE.PlaneGeometry(planeSize, planeSize);
    
    const defaultPlaneNormal = new THREE.Vector3(0, 0, 1); // PlaneGeometry default normal
    const quaternion = new THREE.Quaternion().setFromUnitVectors(defaultPlaneNormal, planeNormalNormalized);

    return { geometry, position: finalPlaneCenterPosition, quaternion };

  }, [showMillerPlane, millerH, millerK, millerL, repeats, unitCellSceneOrigin]); // unitCellSceneOrigin is now a prop

  return (
    <>
      {atoms.map((atom, index) => (
        <mesh key={`atom-${index}-${atom.x}-${atom.y}-${atom.z}`} position={[atom.x, atom.y, atom.z]}>
          <sphereGeometry args={[atom.radius, 32, 32]} />
          <meshStandardMaterial color={atom.color} metalness={0.3} roughness={0.6} />
        </mesh>
      ))}

      {unitCellEdges && (
        <lineSegments geometry={unitCellEdges}>
          <lineBasicMaterial color={0xffffff} linewidth={1.5} />
        </lineSegments>
      )}

      {millerPlaneData && (
        <mesh
          geometry={millerPlaneData.geometry}
          position={millerPlaneData.position}
          quaternion={millerPlaneData.quaternion}
        >
          <meshStandardMaterial 
            color={0x38A169} // A slightly different green
            side={THREE.DoubleSide} 
            transparent 
            opacity={0.55}
            metalness={0.1}
            roughness={0.7}
          />
        </mesh>
      )}
    </>
  );
}

export default function CrystalVisualizationPage() {
  const [crystalStructure, setCrystalStructure] = useState("SC");
  const [repeats, setRepeats] = useState(1); // Default to 1 from screenshot
  const [radiusPrimary, setRadiusPrimary] = useState(0.071); // From screenshot
  const [radiusSecondary, setRadiusSecondary] = useState(0.215); // From screenshot
  const [shiftX, setShiftX] = useState(-0.06); // From screenshot
  const [shiftY, setShiftY] = useState(0.00);
  const [shiftZ, setShiftZ] = useState(0.00);
  const [showUnitCell, setShowUnitCell] = useState(true);
  const [showMillerPlane, setShowMillerPlane] = useState(true); // From screenshot
  const [millerH, setMillerH] = useState("1");
  const [millerK, setMillerK] = useState("1");
  const [millerL, setMillerL] = useState("1");

  // Force re-mount of Canvas for camera re-calculation if repeats change significantly
  const canvasKey = useMemo(() => `${repeats}-${crystalStructure}`, [repeats, crystalStructure]); 

  const unitCellSceneOrigin = useMemo(() => {
    const supercellCenteringOffset = (repeats > 1) ? (repeats * latticeConstant - latticeConstant) / 2 : 0;
    return new THREE.Vector3(
        -supercellCenteringOffset,
        -supercellCenteringOffset,
        -supercellCenteringOffset
    );
  }, [repeats]);

  return (
    <div style={{ display: 'flex', height: 'calc(100vh - 40px)', margin: '20px', gap: '20px' }}>
      <Card className="w-[350px] min-w-[300px] shadow-lg">
        <CardHeader>
          <CardTitle>Controls</CardTitle>
        </CardHeader>
        <CardContent className="space-y-6 overflow-y-auto" style={{maxHeight: 'calc(100vh - 120px)'}}> {/* Adjusted maxHeight */}
          <div>
            <Label htmlFor="crystal-structure">Crystal Structure</Label>
            <Select value={crystalStructure} onValueChange={setCrystalStructure}>
              <SelectTrigger id="crystal-structure">
                <SelectValue placeholder="Select structure" />
              </SelectTrigger>
              <SelectContent>
                <SelectItem value="SC">Simple Cubic (SC)</SelectItem>
                <SelectItem value="BCC">Body-Centered Cubic (BCC)</SelectItem>
                <SelectItem value="FCC">Face-Centered Cubic (FCC)</SelectItem>
              </SelectContent>
            </Select>
          </div>

          <div>
            <Label htmlFor="repeats" className="flex justify-between">
              <span>Repeats</span> <span>({repeats})</span>
            </Label>
            <Slider
              id="repeats"
              min={1} max={5} step={1}
              value={[repeats]}
              onValueChange={(val) => setRepeats(val[0])}
            />
          </div>

          <div>
            <Label htmlFor="radius-primary" className="flex justify-between">
                <span>Radius (Primary)</span><span>({radiusPrimary.toFixed(3)})</span>
            </Label>
            <Slider
              id="radius-primary"
              min={0.01} max={0.75} step={0.001} // Min radius smaller
              value={[radiusPrimary]}
              onValueChange={(val) => setRadiusPrimary(val[0])}
            />
          </div>
          
          <div>
            <Label htmlFor="radius-secondary" className="flex justify-between">
                <span>Radius (Secondary)</span><span>({radiusSecondary.toFixed(3)})</span>
            </Label>
            <Slider
              id="radius-secondary"
              min={0.01} max={0.75} step={0.001}
              value={[radiusSecondary]}
              onValueChange={(val) => setRadiusSecondary(val[0])}
              disabled // Keep disabled unless specific use case for it arises
            />
          </div>

          <div className="space-y-1">
            <Label>Cell Shifts (Fractional)</Label>
            <div>
              <Label htmlFor="shift-x" className="text-xs flex justify-between">
                <span>Shift X</span><span>({shiftX.toFixed(2)})</span>
              </Label>
              <Slider id="shift-x" min={-0.5} max={0.5} step={0.01} value={[shiftX]} onValueChange={(val) => setShiftX(val[0])} />
            </div>
            <div>
              <Label htmlFor="shift-y" className="text-xs flex justify-between">
                <span>Shift Y</span><span>({shiftY.toFixed(2)})</span>
              </Label>
              <Slider id="shift-y" min={-0.5} max={0.5} step={0.01} value={[shiftY]} onValueChange={(val) => setShiftY(val[0])} />
            </div>
            <div>
              <Label htmlFor="shift-z" className="text-xs flex justify-between">
                <span>Shift Z</span><span>({shiftZ.toFixed(2)})</span>
              </Label>
              <Slider id="shift-z" min={-0.5} max={0.5} step={0.01} value={[shiftZ]} onValueChange={(val) => setShiftZ(val[0])} />
            </div>
          </div>

          <div className="flex items-center space-x-2">
            <Checkbox id="show-unit-cell" checked={showUnitCell} onCheckedChange={(checked) => setShowUnitCell(Boolean(checked))} />
            <Label htmlFor="show-unit-cell">Show Unit Cell</Label>
          </div>
          
          <div className="space-y-2 pt-2 border-t">
            <Label className="text-base font-semibold">Miller Indices (hkl)</Label>
             <div className="flex items-center space-x-2">
                <Checkbox id="show-miller-plane" checked={showMillerPlane} onCheckedChange={(checked) => setShowMillerPlane(Boolean(checked))} />
                <Label htmlFor="show-miller-plane">Show Miller Plane(s)</Label>
            </div>
            <div className="grid grid-cols-3 gap-2 items-end">
              <div>
                <Label htmlFor="miller-h" className="text-xs">h</Label>
                <Input id="miller-h" type="number" value={millerH} onChange={(e) => setMillerH(e.target.value)} placeholder="h"/>
              </div>
              <div>
                <Label htmlFor="miller-k" className="text-xs">k</Label>
                <Input id="miller-k" type="number" value={millerK} onChange={(e) => setMillerK(e.target.value)} placeholder="k"/>
              </div>
              <div>
                <Label htmlFor="miller-l" className="text-xs">l</Label>
                <Input id="miller-l" type="number" value={millerL} onChange={(e) => setMillerL(e.target.value)} placeholder="l"/>
              </div>
            </div>
          </div>

        </CardContent>
      </Card>

      <div style={{ flex: 1, border: '1px solid hsl(var(--border))', borderRadius: 'var(--radius)', overflow: 'hidden', background: '#2D3748' /* A dark bg */ }}>
        <Canvas
          key={canvasKey}
          camera={{ 
            position: [latticeConstant * repeats * 1.6, latticeConstant * repeats * 1.6, latticeConstant * repeats * 1.6], 
            fov: 50,
            near: 0.1,
            far: latticeConstant * repeats * 10
          }}
          shadows
        >
          <Suspense fallback={<Html center><p style={{color: 'white'}}>Loading 3D Scene...</p></Html>}>
            <ambientLight intensity={Math.PI / 2 * 0.6} /> {/* Adjusted intensity based on common three examples */}
            <directionalLight 
              position={[latticeConstant * repeats * 1.5, latticeConstant * repeats * 2.5, latticeConstant * repeats * 2]} 
              intensity={Math.PI * 0.8}
              castShadow 
              shadow-mapSize-width={1024}
              shadow-mapSize-height={1024}
            />
            <CrystalScene
              crystalStructure={crystalStructure}
              repeats={repeats}
              radiusPrimary={radiusPrimary}
              radiusSecondary={radiusSecondary}
              shiftX={shiftX}
              shiftY={shiftY}
              shiftZ={shiftZ}
              showUnitCell={showUnitCell}
              showMillerPlane={showMillerPlane}
              millerH={millerH}
              millerK={millerK}
              millerL={millerL}
              unitCellSceneOrigin={unitCellSceneOrigin} // Pass as prop
            />
            <OrbitControls 
              enableDamping 
              dampingFactor={0.05}
              minDistance={latticeConstant * 0.5}
              maxDistance={latticeConstant * repeats * 6} // Increased max distance slightly
            />
            <gridHelper 
                args={[latticeConstant * repeats * 2.5, repeats * 10, '#555555', '#777777']} 
                position={[0, unitCellSceneOrigin.y - latticeConstant * 0.01, 0]} // Use lifted state
            />
          </Suspense>
        </Canvas>
      </div>
    </div>
  );
}