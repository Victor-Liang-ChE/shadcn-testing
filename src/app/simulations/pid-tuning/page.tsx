"use client";
import React, { useState, useEffect, useCallback, useMemo, useRef } from 'react';

// ECharts imports
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import type { EChartsOption, SeriesOption, EChartsType } from 'echarts';
import { LineChart } from 'echarts/charts';
import { TitleComponent, TooltipComponent, GridComponent, LegendComponent, ToolboxComponent } from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';

echarts.use([TitleComponent, TooltipComponent, GridComponent, LegendComponent, ToolboxComponent, LineChart, CanvasRenderer]);

// Shadcn UI imports
import { Card, CardContent } from "@/components/ui/card";
import { Label } from "@/components/ui/label";
import { Slider } from "@/components/ui/slider";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Skeleton } from "@/components/ui/skeleton";
import { Input } from "@/components/ui/input";

// --- RK45 Solver ---
// (RK45 implementation assumed correct and unchanged)
function solveODERK45(
    dydt: (t: number, y: number[]) => number[],
    y0: number[],
    tSpan: [number, number],
    tol: number = 1e-5,
    h_initial: number = 0.01,
    maxSteps: number = 5000
): { t: number[]; y: number[][]; status: string } {
    // --- Implementation from previous steps ---
    const [t0, tf] = tSpan; let t = t0; let y = y0.slice(); const t_arr = [t]; const y_arr = [y.slice()];
    let h = h_initial; const safety = 0.84; const minFactor = 0.2; const maxFactor = 5.0;
    const c2=1/5, c3=3/10, c4=4/5, c5=8/9, c6=1, c7=1; const a21=1/5; const a31=3/40, a32=9/40;
    const a41=44/45, a42=-56/15, a43=32/9; const a51=19372/6561, a52=-25360/2187, a53=64448/6561, a54=-212/729;
    const a61=9017/3168, a62=-355/33, a63=46732/5247, a64=49/176, a65=-5103/18656;
    const a71=35/384, a72=0, a73=500/1113, a74=125/192, a75=-2187/6784, a76=11/84;
    const b1=35/384, b2=0, b3=500/1113, b4=125/192, b5=-2187/6784, b6=11/84, b7=0; // 5th order (y5th)
    const b1s=5179/57600, b2s=0, b3s=7571/16695, b4s=393/640, b5s=-92097/339200, b6s=187/2100, b7s=1/40; // 4th order (y4th)
    let stepCount = 0; let status = "Complete";
    while (t < tf && stepCount < maxSteps) {
        if (t + h > tf) h = tf - t; if (h <= 1e-15) { status = "Warning: Step size too small."; break; }
        const k1 = dydt(t, y);
        const y2 = y.map((yi: number, i: number) => yi + h * a21 * k1[i]);
        const k2 = dydt(t + c2 * h, y2);
        const y3 = y.map((yi: number, i: number) => yi + h * (a31 * k1[i] + a32 * k2[i]));
        const k3 = dydt(t + c3 * h, y3);
        const y4 = y.map((yi: number, i: number) => yi + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]));
        const k4 = dydt(t + c4 * h, y4);
        const y5 = y.map((yi: number, i: number) => yi + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]));
        const k5 = dydt(t + c5 * h, y5);
        const y6 = y.map((yi: number, i: number) => yi + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]));
        const k6 = dydt(t + c6 * h, y6);
        const y7 = y.map((yi: number, i: number) => yi + h * (a71*k1[i] + a73*k3[i] + a74*k4[i] + a75*k5[i] + a76*k6[i]));
        const k7 = dydt(t+h, y7);

        const y5th = y.map((yi: number, i: number) => yi + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]));
        const y4th = y.map((yi: number, i: number) => yi + h * (b1s * k1[i] + b3s * k3[i] + b4s * k4[i] + b5s * k5[i] + b6s * k6[i] + b7s * k7[i]));

        let error_norm = 0;
        for (let i = 0; i < y.length; i++) {
            const scale = tol + Math.abs(y[i]) * tol;
            const error_i = Math.abs(y5th[i] - y4th[i]);
            error_norm += (error_i / scale) ** 2;
        }
        const error_ratio = Math.sqrt(error_norm / y.length);

        if (error_ratio <= 1 || h < 1e-14) {
            t += h; y = y5th.slice(); t_arr.push(t); y_arr.push(y.slice());
        }
        let factor = safety * Math.pow(Math.max(error_ratio, 1e-10), -0.2);
        factor = Math.min(maxFactor, Math.max(minFactor, factor));
        h = h * factor;
        stepCount++;
    }
    if (stepCount >= maxSteps) status = "Warning: Max steps reached."; else if (t < tf && status === "Complete") status = "Warning: Solver stopped early.";
    const numStates = y0.length; const numTimePoints = t_arr.length; const y_transposed: number[][] = Array.from({ length: numStates }, () => Array(numTimePoints).fill(0));
    for (let i = 0; i < numTimePoints; i++) { for (let j = 0; j < numStates; j++) { if (y_arr[i] && y_arr[i][j] !== undefined) { y_transposed[j][i] = y_arr[i][j]; } } }
    return { t: t_arr, y: y_transposed, status };
}

// --- Type Definitions ---
interface PidParams { Kp: number; Ti: number; Td: number; }
type ProcessParameterValues = { K?: number; T?: number; T1?: number; T2?: number; T3?: number; theta?: number; zeta?: number; beta?: number; };
interface ModelParameter { id: keyof ProcessParameterValues; label: string; symbol: string; unit?: string; default: number; min: number; max: number; step: number; }
type ProcessSimulationResult = { t: number[]; y: number[][]; status: string };
interface ProcessModel { id: string; label: string; descriptionHtml: string; parameters: ModelParameter[]; applicableRuleIds: string[]; getOdeSystem: (params: ProcessParameterValues, setpoint: number, pid: PidParams, alpha: number) => { odeFunc: (t: number, y: number[]) => number[], y0: number[], C_proc: number[][] }; }
interface ImcRule { id: string; label: string; calculate: (params: ProcessParameterValues & { lambda: number }) => PidParams; }


// --- Helper: State Space PID Controller ---
const getPidDerivatives = ( pid: PidParams, alpha: number, setpoint: number, processOutput: number, controllerStates: { xi: number, xf: number } ): { dxi: number; dxf: number; u: number } => {
    const { Kp, Ti, Td } = pid; const { xi, xf } = controllerStates;
    const Ki = (Ti > 1e-9 && isFinite(Ti)) ? Kp / Ti : 0;
    const Tf = (Td > 1e-9) ? Math.max(alpha * Td, 1e-9) : 1e-9;
    const e = setpoint - processOutput;
    const dxi = e;
    const dxf = (Td > 1e-9) ? (1 / Tf) * (e - xf) : 0;
    const derivativeTerm = (Td > 1e-9) ? (Kp * Td / Tf) * (e - xf) : 0;
    const u = Kp * e + Ki * xi + derivativeTerm;
    return { dxi, dxf, u };
};

// --- Process Model Definitions ---
const processModels: ProcessModel[] = [
    // --- Case A: First Order ---
    {
        id: 'case_a', label: 'A: First Order', descriptionHtml: 'K / (&tau;s + 1)',
        parameters: [ { id: 'K', label: 'Gain', symbol:'K', default: 2.0, min: 0.1, max: 10, step: 0.1 }, { id: 'T', label: 'Time Constant', symbol:'&tau;', default: 5.0, min: 0.1, max: 50, step: 0.1 }, ],
        applicableRuleIds: ['imc_case_a', 'zn_pi'],
        getOdeSystem: (p, setpoint, pid, alpha) => {
            const K = p.K ?? 1; const T = Math.max(p.T ?? 1, 1e-9);
            // FIX: Rename A, B, C to avoid potential scope conflicts
            const A_ss = [[-1 / T]]; const B_ss = [[K / T]]; const C_ss = [[1]];
            const order = 1;
            const ode = (t: number, x: number[]): number[] => {
                const xp = x[0]; const xi = x[1]; const xf = x[2];
                const y = C_ss[0][0] * xp;
                const { dxi, dxf, u } = getPidDerivatives(pid, alpha, setpoint, y, { xi, xf });
                const dxp = A_ss[0][0] * xp + B_ss[0][0] * u;
                return [dxp, dxi, dxf];
            };
            const y0 = Array(order + 2).fill(0);
            return { odeFunc: ode, y0, C_proc: C_ss };
        }
    },
    // --- Case B: Second Order (Overdamped) ---
     {
        id: 'case_b', label: 'B: 2nd Order (Overdamped)', descriptionHtml: 'K / [(<span style="font-style: italic">&tau;<sub>1</sub>s + 1)(&tau;<sub>2</sub>s + 1)</span>]',
        parameters: [ { id: 'K', label: 'Gain', symbol:'K', default: 2.0, min: 0.1, max: 10, step: 0.1 }, { id: 'T1', label: 'Time Constant 1', symbol:'&tau;<sub>1</sub>', default: 5.0, min: 0.1, max: 50, step: 0.1 }, { id: 'T2', label: 'Time Constant 2', symbol:'&tau;<sub>2</sub>', default: 2.0, min: 0.1, max: 50, step: 0.1 }, ],
        applicableRuleIds: ['imc_case_b'],
        getOdeSystem: (p, setpoint, pid, alpha) => {
            const K = p.K ?? 1; const T1 = Math.max(p.T1 ?? 1, 1e-9); const T2 = Math.max(p.T2 ?? 1, 1e-9);
            // FIX: Rename A, B, C
            const A_ss = [ [0, 1], [-1/(T1*T2), -(T1+T2)/(T1*T2)] ]; const B_ss = [ [0], [K/(T1*T2)] ]; const C_ss = [[1, 0]];
            const order = 2;
            const ode = (t: number, x: number[]): number[] => {
                const xp = x.slice(0, order); const xi = x[order]; const xf = x[order + 1];
                const y = C_ss[0][0] * xp[0] + C_ss[0][1] * xp[1];
                const { dxi, dxf, u } = getPidDerivatives(pid, alpha, setpoint, y, { xi, xf });
                const dxp1 = A_ss[0][0]*xp[0] + A_ss[0][1]*xp[1] + B_ss[0][0]*u;
                const dxp2 = A_ss[1][0]*xp[0] + A_ss[1][1]*xp[1] + B_ss[1][0]*u;
                return [dxp1, dxp2, dxi, dxf];
            };
            const y0 = Array(order + 2).fill(0);
            return { odeFunc: ode, y0, C_proc: C_ss };
        }
    },
    // --- Case C: Second Order (Underdamped/Critically Damped) ---
    {
        id: 'case_c', label: 'C: 2nd Order (Underdamped)', descriptionHtml: 'K / (<span style="font-style: italic">&tau;<sup>2</sup>s<sup>2</sup> + 2&zeta;&tau;s + 1</span>)',
        parameters: [ { id: 'K', label: 'Gain', symbol:'K', default: 2.0, min: 0.1, max: 10, step: 0.1 }, { id: 'T', label: 'Natural Period', symbol:'&tau;', default: 5.0, min: 0.1, max: 50, step: 0.1 }, { id: 'zeta', label: 'Damping Factor', symbol:'&zeta;', default: 0.5, min: 0.1, max: 5, step: 0.05 }, ],
        applicableRuleIds: ['imc_case_c'],
        getOdeSystem: (p, setpoint, pid, alpha) => {
            const K = p.K ?? 1; const T = Math.max(p.T ?? 1, 1e-9); const zeta = p.zeta ?? 0.5;
            const T_sq = T*T;
            // FIX: Rename A, B, C
            const A_ss = [ [0, 1], [-1/T_sq, -2*zeta/T] ]; const B_ss = [ [0], [K/T_sq] ]; const C_ss = [[1, 0]];
            const order = 2;
            const ode = (t: number, x: number[]): number[] => {
                const xp = x.slice(0, order); const xi = x[order]; const xf = x[order + 1];
                const y = C_ss[0][0] * xp[0] + C_ss[0][1] * xp[1];
                const { dxi, dxf, u } = getPidDerivatives(pid, alpha, setpoint, y, { xi, xf });
                const dxp1 = A_ss[0][0]*xp[0] + A_ss[0][1]*xp[1] + B_ss[0][0]*u;
                const dxp2 = A_ss[1][0]*xp[0] + A_ss[1][1]*xp[1] + B_ss[1][0]*u;
                return [dxp1, dxp2, dxi, dxf];
            };
            const y0 = Array(order + 2).fill(0);
            return { odeFunc: ode, y0, C_proc: C_ss };
        }
    },
    // --- Case E: Pure Integrator ---
    {
        id: 'case_e', label: 'E: Pure Integrator', descriptionHtml: 'K / s',
        parameters: [ { id: 'K', label: 'Gain', symbol:'K', default: 0.5, min: 0.01, max: 10, step: 0.01 }, ],
        applicableRuleIds: ['imc_case_e'],
        getOdeSystem: (p, setpoint, pid, alpha) => {
            const K = p.K ?? 1;
            // FIX: Rename A, B, C
            const A_ss = [[0]]; const B_ss = [[K]]; const C_ss = [[1]];
            const order = 1;
            const ode = (t: number, x: number[]): number[] => {
                const xp = x[0]; const xi = x[1]; const xf = x[2];
                const y = C_ss[0][0] * xp;
                const { dxi, dxf, u } = getPidDerivatives(pid, alpha, setpoint, y, { xi, xf });
                const dxp = B_ss[0][0] * u;
                return [dxp, dxi, dxf];
            };
            const y0 = Array(order + 2).fill(0);
            return { odeFunc: ode, y0, C_proc: C_ss };
        }
    },
    // --- Case F: Integrator + First Order Lag ---
    {
        id: 'case_f', label: 'F: Integrator + Lag', descriptionHtml: 'K / [s(&tau;s + 1)]',
        parameters: [ { id: 'K', label: 'Gain', symbol:'K', default: 0.5, min: 0.01, max: 10, step: 0.01 }, { id: 'T', label: 'Lag Time Constant', symbol:'&tau;', default: 5.0, min: 0.1, max: 50, step: 0.1 }, ],
        applicableRuleIds: ['imc_case_f'],
        getOdeSystem: (p, setpoint, pid, alpha) => {
            const K = p.K ?? 1; const T = Math.max(p.T ?? 1, 1e-9);
            // FIX: Rename A, B, C
            const A_ss = [ [0, 1], [0, -1/T] ]; const B_ss = [ [0], [K/T] ]; const C_ss = [[1, 0]];
            const order = 2;
             const ode = (t: number, x: number[]): number[] => {
                const xp = x.slice(0, order); const xi = x[order]; const xf = x[order + 1];
                const y = C_ss[0][0] * xp[0] + C_ss[0][1] * xp[1];
                const { dxi, dxf, u } = getPidDerivatives(pid, alpha, setpoint, y, { xi, xf });
                const dxp1 = A_ss[0][0]*xp[0] + A_ss[0][1]*xp[1] + B_ss[0][0]*u;
                const dxp2 = A_ss[1][0]*xp[0] + A_ss[1][1]*xp[1] + B_ss[1][0]*u;
                return [dxp1, dxp2, dxi, dxf];
            };
            const y0 = Array(order + 2).fill(0);
            return { odeFunc: ode, y0, C_proc: C_ss };
        }
    },
    // --- Case G/H: FOPDT (First Order Plus Dead Time) ---
    {
        id: 'case_gh', label: 'G/H: FOPDT', descriptionHtml: 'K e<sup>-&theta;s</sup> / (&tau;s + 1)',
        parameters: [ { id: 'K', label: 'Gain', symbol:'K', default: 2.0, min: 0.1, max: 10, step: 0.1 }, { id: 'T', label: 'Time Constant', symbol:'&tau;', default: 5.0, min: 0.1, max: 50, step: 0.1 }, { id: 'theta', label: 'Dead Time', symbol:'&theta;', default: 1.0, min: 0, max: 10, step: 0.05 }, ],
        applicableRuleIds: ['imc_case_g_pi', 'imc_case_h_pid', 'zn_pi'],
        getOdeSystem: (p, setpoint, pid, alpha) => {
            const K = p.K ?? 1; const T = Math.max(p.T ?? 1, 1e-9); const th = p.theta ?? 0;
            const usePade = th > 1e-9;

            // FIX: Rename A, B, C and use consistent definitions
            const A_pade = [[0, 1], [-2/(T*th), -(T+th/2)*2/(T*th)]];
            const B_pade = [[0], [2*K/(T*th)]]; // K included in B
            const C_pade = [[1, -th/2]];      // C modified because K in B? No, C should reflect Padé approx output eqn. y = x1 - th/2 * x1_dot. If x1'=x2, y=x1-th/2*x2.
                                             // Let's revert to previous structure that maybe worked?
            // Previous possibly working structure:
            // const A_pade = [[0, 1], [-2/(T*th), -2*(T+th/2)/(T*th)]];
            // const B_pade = [[0], [1]]; // K NOT in B
            // const C_pade = [[K, -K*th/2]]; // K handled by C

            const A_no_pade = [[-1 / T]]; const B_no_pade = [[K / T]]; const C_no_pade = [[1]]; // K in B here too

            // Let's try K in B for both consistently
            const A_ss = usePade ? A_pade : A_no_pade;
            const B_ss = usePade ? [[0], [2*K/(T*th)]] : B_no_pade; // Use K in B for Padé too
            const C_ss = usePade ? [[1, -th/2]] : C_no_pade; // Output C matches state structure

            const order = usePade ? 2 : 1;

            const ode = (t: number, x: number[]): number[] => {
                const xp = x.slice(0, order); const xi = x[order]; const xf = x[order + 1];
                let y = 0;
                // Calculate process output y = C_ss @ xp
                for (let i = 0; i < C_ss[0].length; i++) {
                     if (i < xp.length && typeof xp[i] === 'number') { y += C_ss[0][i] * xp[i]; }
                     else { console.warn(`State xp[${i}] mismatch in FOPDT output calc.`); y=NaN; break; }
                }
                if (isNaN(y)) return Array(order+2).fill(NaN); // Bail if output calc failed

                const { dxi, dxf, u } = getPidDerivatives(pid, alpha, setpoint, y, { xi, xf });
                const dxp: number[] = Array(order).fill(0);
                // Calculate process state derivatives: dxp = A_ss @ xp + B_ss @ u
                for (let i = 0; i < order; i++) {
                    for (let j = 0; j < A_ss[i].length; j++) {
                        if (j < xp.length && typeof xp[j] === 'number') { dxp[i] += A_ss[i][j] * xp[j]; }
                        else { console.warn(`State xp[${j}] mismatch in FOPDT state derivative calc.`); dxp[i]=NaN; break; }
                    }
                    if(isNaN(dxp[i])) break; // Stop row calculation if error
                    if (B_ss[i] && typeof B_ss[i][0] === 'number') { dxp[i] += B_ss[i][0] * u; }
                    else { console.warn(`Input matrix B_ss[${i}][0] invalid in FOPDT.`); dxp[i]=NaN; }
                }
                 if (dxp.some(isNaN)) return Array(order+2).fill(NaN); // Bail if state calc failed

                return [...dxp, dxi, dxf];
            };
            const y0 = Array(order + 2).fill(0);
            return { odeFunc: ode, y0, C_proc: C_ss };
        }
    },
     // --- Case M: Integrator + Delay ---
    {
        id: 'case_m', label: 'M: Integrator + Delay', descriptionHtml: 'K e<sup>-&theta;s</sup> / s',
        parameters: [ { id: 'K', label: 'Gain', symbol:'K', default: 0.5, min: 0.01, max: 10, step: 0.01 }, { id: 'theta', label: 'Dead Time', symbol:'&theta;', default: 1.0, min: 0, max: 10, step: 0.05 }, ],
        applicableRuleIds: ['imc_case_m_pi', 'imc_case_n_pid'],
        getOdeSystem: (p, setpoint, pid, alpha) => {
            const K = p.K ?? 1; const th = Math.max(p.theta ?? 0, 1e-9);
            // FIX: Rename A, B, C
            // State-space from previous derivation:
            const A_ss = [ [-2/th, 1], [0, 0] ]; const B_ss = [ [-K], [2*K/th] ]; const C_ss = [[1, 0]];
            const order = 2;
            const ode = (t: number, x: number[]): number[] => {
                const xp = x.slice(0, order); const xi = x[order]; const xf = x[order + 1];
                const y = C_ss[0][0] * xp[0] + C_ss[0][1] * xp[1];
                const { dxi, dxf, u } = getPidDerivatives(pid, alpha, setpoint, y, { xi, xf });
                const dxp1 = A_ss[0][0]*xp[0] + A_ss[0][1]*xp[1] + B_ss[0][0]*u;
                const dxp2 = A_ss[1][0]*xp[0] + A_ss[1][1]*xp[1] + B_ss[1][0]*u;
                return [dxp1, dxp2, dxi, dxf];
            };
            const y0 = Array(order + 2).fill(0);
            return { odeFunc: ode, y0, C_proc: C_ss };
        }
    },
];

// --- IMC Tuning Rules (Matching Table 12.1 Cases) ---
const imcRules: ImcRule[] = [
    { id: 'imc_case_a', label: 'Case A: PI (IMC)', calculate: ({ K=1, T=1, lambda }) => { const safe_lambda = Math.max(lambda, 1e-9); const Kp = T / (K * safe_lambda); const Ti = T; return { Kp, Ti: Math.max(Ti, 1e-9), Td: 0 }; } },
    { id: 'imc_case_b', label: 'Case B: PI (IMC)', calculate: ({ K=1, T1=1, T2=1, lambda }) => { const safe_lambda = Math.max(lambda, 1e-9); const Kp = (T1 + T2) / (K * safe_lambda); const Ti = T1 + T2; return { Kp, Ti: Math.max(Ti, 1e-9), Td: 0 }; } },
    { id: 'imc_case_c', label: 'Case C: PI (IMC)', calculate: ({ K=1, T=1, zeta=1, lambda }) => { const safe_lambda = Math.max(lambda, 1e-9); const Kp = (2 * zeta * T) / (K * safe_lambda); const Ti = 2 * zeta * T; return { Kp, Ti: Math.max(Ti, 1e-9), Td: 0 }; } },
    { id: 'imc_case_e', label: 'Case E: PI (IMC)', calculate: ({ K=1, lambda }) => { const safe_lambda = Math.max(lambda, 1e-9); const Kp = 2 / (K * safe_lambda); const Ti = 2 * safe_lambda; return { Kp, Ti: Math.max(Ti, 1e-9), Td: 0 }; } },
    { id: 'imc_case_f', label: 'Case F: PID (IMC)', calculate: ({ K=1, T=1, lambda }) => {
        const safe_lambda = Math.max(lambda, 1e-9);
        const T_eff = Math.max(T, 1e-9); // Use effective T
        // FIX: Declare Kp, Ti, Td with const
        const Kp = (2*T_eff + safe_lambda) / (K * (safe_lambda**2));
        const Ti = 2*T_eff + safe_lambda;
        const Td = (2*T_eff*safe_lambda) / (2*T_eff + safe_lambda);
        return { Kp, Ti: Math.max(Ti, 1e-9), Td: Math.max(Td, 0) };
    }},
    { id: 'imc_case_g_pi', label: 'Case G: PI (IMC)', calculate: ({ K=1, T=1, theta=0, lambda }) => { const safe_lambda = Math.max(lambda, theta, 1e-9); const T_eff = Math.max(T, 1e-9); const Kp = T_eff / (K * (safe_lambda + theta)); const Ti = T_eff; return { Kp, Ti: Math.max(Ti, 1e-9), Td: 0 }; } },
    { id: 'imc_case_h_pid', label: 'Case H: PID (IMC)', calculate: ({ K=1, T=1, theta=0, lambda }) => { const T_eff = Math.max(T, 1e-9); const th_eff = Math.max(theta, 1e-9); const safe_lambda = Math.max(lambda, 0.5*th_eff, 1e-9); const Kp = (T_eff + th_eff / 2) / (K * (safe_lambda + th_eff / 2)); const Ti = T_eff + th_eff / 2; const Td = (T_eff * th_eff) / (2 * T_eff + th_eff); return { Kp, Ti: Math.max(Ti, 1e-9), Td: Math.max(Td, 0) }; } },
    { id: 'imc_case_m_pi', label: 'Case M: PI (IMC)', calculate: ({ K=1, theta=0, lambda }) => { const th_eff = Math.max(theta, 1e-9); const safe_lambda = Math.max(lambda, th_eff, 1e-9); const Kp = 2 / (K * (2*safe_lambda + th_eff)); const Ti = 2*safe_lambda + th_eff; return { Kp, Ti: Math.max(Ti, 1e-9), Td: 0 }; } },
    { id: 'imc_case_n_pid', label: 'Case N: PID (IMC)', calculate: ({ K=1, theta=0, lambda }) => { const th_eff = Math.max(theta, 1e-9); const safe_lambda = Math.max(lambda, 0.5*th_eff, 1e-9); const denom = K * (safe_lambda + th_eff/2)**2; const Kp = (2 * safe_lambda + th_eff) / denom; const Ti = safe_lambda + th_eff / 2; const Td = (safe_lambda * th_eff) / (2*safe_lambda + th_eff); return { Kp, Ti: Math.max(Ti, 1e-9), Td: Math.max(Td, 0) }; } },
    { id: 'zn_pi', label: 'Ziegler-Nichols PI', calculate: (params) => { const { K = 1, T = 1, theta = 0.1 } = params; const T_eff = Math.max(T, 1e-9); const safe_theta = Math.max(theta, 1e-9); const Kp = 0.9 * T_eff / (K * safe_theta); const Ti = safe_theta / 0.3; return { Kp, Ti: Math.max(Ti, 1e-9), Td: 0 }; } },
];


// --- PID Tuning Page Component ---
export default function PidTuningPage() {
    const alpha = 0.1; const setpoint = 1.0;
    const [tFinal, setTFinal] = useState<number>(30.0);
    const [selectedModelId, setSelectedModelId] = useState<string>(processModels[0].id);
    const [modelParams, setModelParams] = useState<ProcessParameterValues>({});
    const [applicableRules, setApplicableRules] = useState<ImcRule[]>([]);
    const [selectedRuleId, setSelectedRuleId] = useState<string>('');
    const [lambda, setLambda] = useState<number>(1.5);
    const [pidParams, setPidParams] = useState<PidParams | null>(null);
    const [simulationResult, setSimulationResult] = useState<ProcessSimulationResult | null>(null);
    const [simulationStatus, setSimulationStatus] = useState<string>("Ready");
    const [isLoading, setIsLoading] = useState<boolean>(true);
    const [echartsOptions, setEchartsOptions] = useState<EChartsOption>({});
    const [currentCProc, setCurrentCProc] = useState<number[][]>([[1]]);
    const echartsRef = useRef<ReactECharts | null>(null);
    const selectedModel = useMemo(() => processModels.find(m => m.id === selectedModelId) || processModels[0], [selectedModelId]);
    const selectedRule = useMemo(() => applicableRules.find(rule => rule.id === selectedRuleId), [selectedRuleId, applicableRules]);

    // --- Effects ---
    useEffect(() => {
        const initialParams: ProcessParameterValues = {};
        selectedModel.parameters.forEach(p => { initialParams[p.id] = modelParams[p.id] !== undefined ? modelParams[p.id] : p.default; });
        setModelParams(initialParams);
        const rules = imcRules.filter(rule => selectedModel.applicableRuleIds.includes(rule.id));
        setApplicableRules(rules);
        const firstApplicableRuleId = rules[0]?.id ?? '';
        setSelectedRuleId(firstApplicableRuleId);
        setSimulationResult(null); setPidParams(null); setSimulationStatus("Ready"); setIsLoading(true);
    }, [selectedModelId]);

    // --- Callbacks ---
    const handleModelParamChange = useCallback((paramId: keyof ProcessParameterValues, value: number | string) => {
        const paramDef = selectedModel.parameters.find(p => p.id === paramId);
        let numericValue = typeof value === 'string' ? parseFloat(value) : value;
        if (isNaN(numericValue)) numericValue = paramDef?.default ?? 0;
        const validatedValue = Math.min(Math.max(numericValue, paramDef?.min ?? -Infinity), paramDef?.max ?? Infinity);
        setModelParams(prev => ({ ...prev, [paramId]: validatedValue }));
    }, [selectedModel.parameters]);

    const formatNumber = (num: number | null | undefined, precision = 3): string => { /* ... (unchanged) ... */ if (num === null || num === undefined || isNaN(num)) return "---"; if (Math.abs(num) === Infinity) return "Inf"; if (Math.abs(num) > 0 && (Math.abs(num) < 10**(-precision) || Math.abs(num) > 10**6)) return num.toPrecision(precision); return num.toFixed(precision); };

    const generateEchartsOptions = useCallback(( simRes: ProcessSimulationResult | null, calculatedPid: PidParams | null, ruleLabel: string | undefined, currentLambda: number, modelLabel: string, C_proc_used: number[][], currentTFinal: number ): EChartsOption => { /* ... (unchanged from previous version) ... */ const subTextParts = []; if (ruleLabel) subTextParts.push(`Rule: ${ruleLabel}`); if (calculatedPid) { subTextParts.push(`&tau;<sub>c</sub>=${currentLambda.toFixed(2)}`); subTextParts.push(`Kp=${formatNumber(calculatedPid.Kp, 2)}`); subTextParts.push(`&tau;<sub>i</sub>=${formatNumber(calculatedPid.Ti, 2)}`); if (calculatedPid.Td > 1e-9) { subTextParts.push(`&tau;<sub>d</sub>=${formatNumber(calculatedPid.Td, 2)}`); } } const subtext = subTextParts.join(', '); if (!simRes || !C_proc_used || C_proc_used.length === 0 || C_proc_used[0].length === 0) { return { title: { text: 'Waiting for simulation data...', subtext: 'Adjust parameters or select rule', left: 'center' }, xAxis: { type: 'value', min:0, max:currentTFinal }, yAxis: { type: 'value' }, series: [], backgroundColor: 'transparent' }; } const numProcessStates = C_proc_used[0].length; if (numProcessStates > simRes.y.length) { console.error("C Matrix columns mismatch simulation states rows!", { C_cols: numProcessStates, y_states_rows: simRes.y.length }); return { title: { text: 'Error: State/Output matrix row mismatch', left: 'center' }, backgroundColor: 'transparent' }; } const output_y = simRes.t.map((_, i) => { let y_val = 0; if (simRes.y && simRes.y[0] && i < simRes.y[0].length) { for (let j = 0; j < numProcessStates; j++) { if (simRes.y[j] && typeof simRes.y[j][i] === 'number' && C_proc_used[0] && typeof C_proc_used[0][j] === 'number') { y_val += C_proc_used[0][j] * simRes.y[j][i]; } else { console.warn(`State y[${j}][${i}] or C[0][${j}] invalid.`); return NaN; } } } else { console.warn(`Time index ${i} out of bounds.`); return NaN; } return y_val; }); const validDataPoints = simRes.t.map((time, index) => [time, output_y[index]]).filter(point => !isNaN(point[1] as number)); const seriesData: SeriesOption[] = [ { name: 'Process Output (y)', type: 'line', smooth: true, symbol: 'none', data: validDataPoints, lineStyle: { color: '#3b82f6', width: 2 }, emphasis: { focus: 'series' }, }, { name: 'Setpoint (r)', type: 'line', smooth: false, symbol: 'none', data: [[0, setpoint], [currentTFinal, setpoint]], lineStyle: { color: '#ef4444', type: 'dashed' }, step: 'end', z: 1 } ]; return { backgroundColor: 'transparent', animation: false, title: { text: `Closed-Loop Step Response: ${modelLabel}`, subtext, left: 'center', textStyle: { fontSize: 16 }, subtextStyle: { fontSize: 11, color: '#6b7280' } }, tooltip: { trigger: 'axis', formatter: (params: any) => { let t = params[0] ? `Time: ${params[0].axisValueLabel}<br/>` : ''; params.forEach((p: any) => { if (p.seriesName === 'Setpoint (r)') { t += `${p.marker}${p.seriesName}: ${setpoint.toFixed(3)}<br/>`; } else if (p.value && p.value.length > 1 && typeof p.value[1] === 'number') { t += `${p.marker}${p.seriesName}: ${p.value[1].toPrecision(4)}<br/>`; } }); return t; } }, legend: { data: ['Process Output (y)', 'Setpoint (r)'], bottom: 10 }, grid: { left: '8%', right: '8%', bottom: '15%', top: '22%', containLabel: false }, toolbox: { feature: { saveAsImage: { name: `pid_tuning_${modelLabel.replace(/[^a-z0-9]/gi, '_')}_${ruleLabel?.replace(/[^a-z0-9]/gi, '_') ?? 'no_rule'}_tc_${currentLambda.toFixed(2)}` } }, show: true, orient: 'horizontal', bottom: 5, right: '5%' }, xAxis: { type: 'value', name: 'Time (s)', nameLocation: 'middle', nameGap: 25, min: 0, max: currentTFinal, axisLabel: { formatter: (v: number) => v.toFixed(1) }, splitLine: { show: false } }, yAxis: { type: 'value', name: 'Output', nameLocation: 'middle', nameGap: 45, axisLabel: { formatter: (v: number) => v.toPrecision(3) }, splitLine: { show: true, lineStyle: { type: 'dashed', color: '#e5e7eb' } } }, series: seriesData }; }, [setpoint]); // Removed args that are passed in, added setpoint dependency

    const runSimulation = useCallback(async () => { /* ... (unchanged simulation logic from previous version, relies on selectedRule, modelParams etc.) ... */ if (!selectedRule) { setIsLoading(false); setSimulationStatus("Select a tuning rule."); const C_proc_prev = (currentCProc && currentCProc.length > 0) ? currentCProc : [[1]]; const emptyOptions = generateEchartsOptions(null, null, undefined, lambda, selectedModel.label, C_proc_prev, tFinal); setEchartsOptions(emptyOptions); return; } const modelParamKeys = selectedModel.parameters.map(p => p.id); const allParamsSet = modelParamKeys.every(key => modelParams[key] !== undefined && modelParams[key] !== null && !isNaN(modelParams[key] as number)); if (!allParamsSet) { setIsLoading(false); setSimulationStatus("Set all model parameters."); return; } setIsLoading(true); setSimulationStatus("Calculating..."); let calculatedPid: PidParams | null = null; try { calculatedPid = selectedRule.calculate({ ...modelParams, lambda }); setPidParams(calculatedPid); } catch (error) { console.error("PID Calculation Error:", error); setSimulationStatus(`Error calculating PID: ${error instanceof Error ? error.message : String(error)}`); setPidParams(null); setSimulationResult(null); const C_proc_prev = (currentCProc && currentCProc.length > 0) ? currentCProc : [[1]]; const errorOptions = generateEchartsOptions(null, null, selectedRule.label, lambda, selectedModel.label, C_proc_prev, tFinal); setEchartsOptions(errorOptions); try { echartsRef.current?.getEchartsInstance().setOption(errorOptions, { notMerge: true }); } catch (e) { console.error("Error setting ECharts option during PID calc error:", e); } setIsLoading(false); return; } if (!calculatedPid) { setSimulationStatus("Error: Failed to obtain PID params."); setIsLoading(false); return; } let sim_results: ProcessSimulationResult | null = null; let finalStatus = "Error: Unknown simulation failure"; let C_proc_used: number[][] = currentCProc; try { const { odeFunc, y0, C_proc } = selectedModel.getOdeSystem(modelParams, setpoint, calculatedPid, alpha); if (!C_proc || C_proc.length === 0 || C_proc[0].length === 0) throw new Error("Invalid C_proc matrix."); C_proc_used = C_proc; setCurrentCProc(C_proc); const t_span: [number, number] = [0, tFinal]; sim_results = await solveODERK45(odeFunc, y0, t_span); setSimulationResult(sim_results); finalStatus = sim_results.status; } catch (error) { console.error("Simulation Error:", error); finalStatus = `Error in simulation: ${error instanceof Error ? error.message : String(error)}`; setSimulationResult(null); setCurrentCProc([[1]]); C_proc_used = [[1]]; } finally { setSimulationStatus(finalStatus); const newOptions = generateEchartsOptions(sim_results, calculatedPid, selectedRule.label, lambda, selectedModel.label, C_proc_used, tFinal); setEchartsOptions(newOptions); try { const instance = echartsRef.current?.getEchartsInstance(); if (instance) instance.setOption(newOptions, { notMerge: true }); else console.warn("ECharts instance not available on ref yet."); } catch (e) { console.error("Error setting ECharts option via ref:", e); } setIsLoading(false); } }, [selectedRule, modelParams, lambda, selectedModel, alpha, setpoint, tFinal, currentCProc, generateEchartsOptions]);

    // Effect to trigger simulation run (debounced)
    useEffect(() => { /* ... (unchanged trigger logic from previous version) ... */ const modelParamKeys = selectedModel.parameters.map(p => p.id); const allParamsSet = modelParamKeys.every(key => modelParams[key] !== undefined && modelParams[key] !== null && !isNaN(modelParams[key] as number)); if (allParamsSet && selectedRuleId && applicableRules.length > 0) { const timer = setTimeout(() => { runSimulation(); }, 300); return () => clearTimeout(timer); } else { setIsLoading(false); } }, [modelParams, selectedRuleId, lambda, selectedModelId, runSimulation, applicableRules, tFinal]);


    // --- Render ---
    const hasValidChartOptions = echartsOptions && echartsOptions.series && Array.isArray(echartsOptions.series);
    return (
        <div className="container mx-auto p-4 md:p-8">
            <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
                {/* Controls Column */}
                <div className="lg:col-span-1 space-y-6">
                     {/* Model Selection */}
                     <Card> <CardContent className="space-y-4 pt-6"> <div> <Label htmlFor="process-model-select" className="mb-2 block font-semibold">Process Model:</Label> <Select value={selectedModelId} onValueChange={(value) => setSelectedModelId(value)} disabled={isLoading}> <SelectTrigger id="process-model-select"><SelectValue placeholder="Select a process model" /></SelectTrigger> <SelectContent> {processModels.map(model => ( <SelectItem key={model.id} value={model.id}> <span className='font-medium'>{model.label}: </span> <span dangerouslySetInnerHTML={{ __html: model.descriptionHtml }} className="ml-2 text-muted-foreground italic" /> </SelectItem> ))} </SelectContent> </Select> </div> </CardContent> </Card>
                     {/* Model Parameters */}
                     <Card> <CardContent className="space-y-4 pt-6"> <Label className="mb-2 block font-semibold">Model Parameters:</Label> {selectedModel.parameters.map(param => ( <div key={param.id} className="space-y-2"> <Label htmlFor={`${param.id}-input`} className="text-sm flex justify-between items-center"> <span dangerouslySetInnerHTML={{ __html: `${param.label} (${param.symbol}):` }} /> <Input id={`${param.id}-input`} type="number" className="h-7 w-24 text-sm px-2 py-1 border rounded" value={modelParams[param.id]?.toString() ?? ''} min={param.min} max={param.max} step={param.step} onChange={(e) => handleModelParamChange(param.id, e.target.value)} disabled={isLoading && simulationStatus === 'Calculating...'} /> </Label> <Slider id={`${param.id}-slider`} min={param.min} max={param.max} step={param.step} value={[modelParams[param.id] ?? param.default]} onValueChange={(v) => handleModelParamChange(param.id, v[0])} disabled={isLoading && simulationStatus === 'Calculating...'} /> </div> ))} </CardContent> </Card>
                     {/* Tuning Rule & Lambda */}
                     <Card> <CardContent className="space-y-4 pt-6"> <div> <Label htmlFor="tuning-rule-select" className="mb-2 block font-semibold">Tuning Rule:</Label> <Select value={selectedRuleId} onValueChange={setSelectedRuleId} disabled={isLoading || applicableRules.length === 0}> <SelectTrigger id="tuning-rule-select"> <SelectValue placeholder={applicableRules.length > 0 ? "Select an applicable rule..." : "No rules applicable..."} /> </SelectTrigger> <SelectContent> {applicableRules.map(rule => ( <SelectItem key={rule.id} value={rule.id}>{rule.label}</SelectItem> ))} </SelectContent> </Select> {applicableRules.length === 0 && !isLoading && <p className="text-xs text-muted-foreground mt-1">No specific IMC rules defined for this model structure yet.</p>} </div> <div className="space-y-2"> <Label htmlFor="lambda-slider" className="text-sm flex justify-between items-center"> <span>IMC Tuning (&tau;<sub>c</sub>):</span> <span className="font-medium">{lambda.toFixed(2)}</span> </Label> <Slider id="lambda-slider" min={0.1} max={Math.max(5.0,(modelParams?.T ?? modelParams?.T1 ?? 1)*1.5)} step={0.05} value={[lambda]} onValueChange={(v)=>setLambda(v[0])} disabled={isLoading && simulationStatus === 'Calculating...'}/> <p className="text-xs text-muted-foreground pt-1">Adjust closed-loop time constant (&tau;<sub>c</sub>).</p> </div> </CardContent> </Card>
                     {/* Calculated PID Parameters */}
                     <Card> <CardContent className="space-y-2 text-sm pt-6 min-h-[100px]"> <Label className="mb-2 block font-semibold">Calculated PID Parameters:</Label> {isLoading && simulationStatus === 'Calculating...' ? ( <> <Skeleton className="h-4 w-20 mb-1"/> <Skeleton className="h-4 w-24 mb-1"/> <Skeleton className="h-4 w-28 mb-1"/> <Skeleton className="h-3 w-32 mt-2"/></> ) : pidParams ? ( <> <p>Kp: <span className="font-semibold">{formatNumber(pidParams.Kp)}</span></p> <p dangerouslySetInnerHTML={{__html:`&tau;<sub>i</sub>: <span class="font-semibold">${formatNumber(pidParams.Ti)}</span>`}}/> {pidParams.Td > 1e-9 && <p dangerouslySetInnerHTML={{__html:`&tau;<sub>d</sub>: <span class="font-semibold">${formatNumber(pidParams.Td)}</span>`}}/>} {pidParams.Td > 1e-9 && <p className="text-xs text-muted-foreground mt-2" dangerouslySetInnerHTML={{__html:`(Filter T<sub>f</sub> &approx; ${formatNumber(alpha*pidParams.Td,2)})`}}/>} </> ) : ( <p className="text-muted-foreground">{simulationStatus.startsWith('Error') ? 'Calculation failed.' : !selectedRuleId ? 'Select tuning rule.' : 'Pending calculation...'}</p> )} </CardContent> </Card>
                     {/* Simulation Time Control */}
                     <Card> <CardContent className="pt-6"> <div className="space-y-2"> <Label htmlFor="tfinal-input" className="text-sm flex justify-between items-center"> <span>Simulation Time (s):</span> <Input id="tfinal-input" type="number" className="h-7 w-24 text-sm px-2 py-1 border rounded" value={tFinal.toString()} min={1} max={500} step={1} onChange={(e) => setTFinal(Math.max(1, parseFloat(e.target.value) || 10))} disabled={isLoading && simulationStatus === 'Calculating...'} /> </Label> <Slider id="tfinal-slider" min={5} max={200} step={5} value={[tFinal]} onValueChange={(v)=>setTFinal(v[0])} disabled={isLoading && simulationStatus === 'Calculating...'}/> </div> </CardContent> </Card>
                </div>
                {/* Plot Column */}
                <div className="lg:col-span-2">
                    <Card> <CardContent className="pt-6"> <div className="text-xs text-muted-foreground mb-2 h-4 text-right pr-2"> {(isLoading && simulationStatus !== 'Calculating...') ? 'Initializing...' : simulationStatus !== 'Complete' && simulationStatus !== 'Ready' ? simulationStatus : ''} </div> <div className="relative h-[500px] md:h-[600px] rounded-md overflow-hidden border"> {isLoading && simulationStatus === 'Calculating...' && ( <div className="absolute inset-0 flex items-center justify-center bg-background/60 dark:bg-background/70 z-10 backdrop-blur-sm" aria-label="Calculating simulation"> <Skeleton className="h-3/4 w-3/4 rounded-md"/> </div> )} {hasValidChartOptions ? ( <ReactECharts ref={echartsRef} echarts={echarts} option={echartsOptions} style={{ height: '100%', width: '100%' }} notMerge={true} lazyUpdate={false} /> ) : ( <div className="absolute inset-0 flex items-center justify-center text-muted-foreground p-4 text-center"> Chart options invalid. </div> )} {!isLoading && !simulationResult && !(simulationStatus.startsWith('Error')) && ( <div className="absolute inset-0 flex items-center justify-center text-muted-foreground p-4 text-center"> { !selectedRuleId ? "Select model and tuning rule to start." : "Adjust parameters to simulate." } </div> )} {simulationStatus.startsWith('Error') && ( <div className="absolute inset-0 flex items-center justify-center text-destructive-foreground bg-destructive/80 p-4 text-center z-20"> {simulationStatus} </div> )} </div> </CardContent> </Card>
                </div>
            </div>
        </div>
    );
}