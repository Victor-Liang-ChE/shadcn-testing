'use client';

import React, { useMemo, Suspense } from 'react';
import * as THREE from 'three';
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Slider } from "@/components/ui/slider";
import { Label } from "@/components/ui/label";
import { Checkbox } from "@/components/ui/checkbox";
import { Input } from "@/components/ui/input";
import { Canvas } from '@react-three/fiber';
import { OrbitControls, Html } from '@react-three/drei';

const latticeConstant = 1.0;

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
  radiusSecondary: number; // Stays in props, though might not be used in CrystalScene
  shiftX: number;
  shiftY: number;
  shiftZ: number;
  showUnitCell: boolean;
  showMillerPlane: boolean;
  millerH: string;
  millerK: string;
  millerL: string;
  // unitCellSceneOrigin removed from props
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
  millerL
}: CrystalSceneProps) {

  const { atoms, unitCellSceneOrigin } = useMemo(() => {
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

    const supercellCenteringOffset = (repeats > 1) ? (repeats * latticeConstant - latticeConstant) / 2 : 0;
    const calculatedUnitCellSceneOrigin = new THREE.Vector3(
        -supercellCenteringOffset,
        -supercellCenteringOffset,
        -supercellCenteringOffset
    );

    for (let rX = 0; rX < repeats; rX++) {
      for (let rY = 0; rY < repeats; rY++) {
        for (let rZ = 0; rZ < repeats; rZ++) {
          const cellCornerX_ideal = rX * latticeConstant;
          const cellCornerY_ideal = rY * latticeConstant;
          const cellCornerZ_ideal = rZ * latticeConstant;
          basePositions.forEach(bp => {
            // Simplify atom coloring: all primary atoms are the same prominent color
            const atomColor = '#FF6B6B'; 

            generatedAtoms.push({
              x: cellCornerX_ideal + (bp.x + shiftX) * latticeConstant - supercellCenteringOffset,
              y: cellCornerY_ideal + (bp.y + shiftY) * latticeConstant - supercellCenteringOffset,
              z: cellCornerZ_ideal + (bp.z + shiftZ) * latticeConstant - supercellCenteringOffset,
              radius: radiusPrimary,
              color: atomColor
            });
          });
        }
      }
    }
    return { atoms: generatedAtoms, unitCellSceneOrigin: calculatedUnitCellSceneOrigin };
  }, [crystalStructure, repeats, radiusPrimary, shiftX, shiftY, shiftZ]);

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
  }, [showUnitCell, unitCellSceneOrigin]);

  const millerPlaneData = useMemo(() => {
    if (!showMillerPlane) return null;
  
    const h = parseInt(millerH), k = parseInt(millerK), l = parseInt(millerL);
    if (isNaN(h) || isNaN(k) || isNaN(l) || (h === 0 && k === 0 && l === 0)) return null;
  
    // 1.  Geometry big enough to cover the super-cell -------------------
    const size      = repeats * latticeConstant * 3;
    const geometry  = new THREE.PlaneGeometry(size, size);
  
    // 2.  Quaternion that fixes the in-plane spin -----------------------
    const n = new THREE.Vector3(h, k, l).normalize();
  
    const almostParallel = (a: THREE.Vector3, b: THREE.Vector3) => Math.abs(a.dot(b)) > 0.999;
    let ref = new THREE.Vector3(0, 0, 1);
    if (almostParallel(n, ref)) ref.set(0, 1, 0);
    if (almostParallel(n, ref)) ref.set(1, 0, 0);
  
    const u = new THREE.Vector3().crossVectors(ref, n).normalize();
    const v = new THREE.Vector3().crossVectors(n, u); // v is already normalized if u and n are ortho-normalized
  
    const quaternion = new THREE.Quaternion()
          .setFromRotationMatrix(new THREE.Matrix4().makeBasis(u, v, n));
  
    // 3.  Position: “centroid” rule (keeps every plane inside the cell) --
    const a = latticeConstant;
    const center = new THREE.Vector3(
        h ? a / h : a * 0.5,
        k ? a / k : a * 0.5,
        l ? a / l : a * 0.5
    );
    const position = unitCellSceneOrigin.clone().add(center);
  
    return { geometry, position, quaternion };
  }, [
    showMillerPlane, millerH, millerK, millerL,
    repeats, unitCellSceneOrigin // latticeConstant is a module-level const, unitCellSceneOrigin depends on repeats
  ]);

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
  const [crystalStructure, setCrystalStructure] = React.useState("FCC");  // Changed to FCC
  const [repeats, setRepeats] = React.useState(1); 
  const [radiusPrimary, setRadiusPrimary] = React.useState(0.15); // Adjusted
  const [radiusSecondary, setRadiusSecondary] = React.useState(0.1); // Adjusted
  const [shiftX, setShiftX] = React.useState(0.00); 
  const [shiftY, setShiftY] = React.useState(0.00);
  const [shiftZ, setShiftZ] = React.useState(0.00);
  const [showUnitCell, setShowUnitCell] = React.useState(true);
  const [showMillerPlane, setShowMillerPlane] = React.useState(true); 
  const [millerH, setMillerH] = React.useState("1"); // Defaulting to 101
  const [millerK, setMillerK] = React.useState("0");
  const [millerL, setMillerL] = React.useState("1");

  // Force re-mount of Canvas for camera re-calculation if repeats change significantly
  const canvasKey = React.useMemo(() => `${repeats}-${crystalStructure}`, [repeats, crystalStructure]); 

  // Calculate gridHelper Y position based on repeats, as unitCellSceneOrigin is internal to CrystalScene
  // unitCellSceneOrigin.y is -supercellCenteringOffset
  const gridHelperYPosition = React.useMemo(() => {
    const supercellCenteringOffset = (repeats > 1) ? (repeats * latticeConstant - latticeConstant) / 2 : 0;
    return -supercellCenteringOffset - latticeConstant * 0.01;
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
              min={0.01} max={0.75} step={0.001} 
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
              disabled 
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
            <ambientLight intensity={Math.PI} /> 
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
              // unitCellSceneOrigin prop removed
            />
            <OrbitControls 
              enableDamping 
              dampingFactor={0.05}
              minDistance={latticeConstant * 0.5}
              maxDistance={latticeConstant * repeats * 6} // Increased max distance slightly
            />
            <gridHelper 
                args={[latticeConstant * repeats * 2.5, repeats * 10, '#555555', '#777777']} 
                position={[0, gridHelperYPosition, 0]} // Use calculated Y position
            />
          </Suspense>
        </Canvas>
      </div>
    </div>
  );
}