'use client'

import { useEffect, useMemo, useRef } from 'react'

type Props = {
  className?: string
  isDark?: boolean
  particleCount?: number
}

type Particle = {
  x: number
  y: number
  vx: number
  vy: number
  r: number
  type: 0 | 1
}

function prefersReducedMotion(): boolean {
  if (typeof window === 'undefined') return false
  return window.matchMedia?.('(prefers-reduced-motion: reduce)')?.matches ?? false
}

export default function MolecularDynamicsThumbnail({
  className,
  isDark = true,
  particleCount = 28,
}: Props) {
  const canvasRef = useRef<HTMLCanvasElement | null>(null)
  const rafRef = useRef<number | null>(null)
  const particlesRef = useRef<Particle[]>([])
  const runningRef = useRef<boolean>(false)

  const palette = useMemo(() => {
    // Match the MD page particle palette: A=blue, B=red.
    return isDark
      ? {
          bg: 'rgba(0,0,0,0)',
          box: 'rgba(255,255,255,0.22)',
          a: 'rgba(96,165,250,0.95)',
          b: 'rgba(248,113,113,0.95)',
        }
      : {
          bg: 'rgba(0,0,0,0)',
          box: 'rgba(17,24,39,0.22)',
          a: 'rgba(59,130,246,0.92)',
          b: 'rgba(239,68,68,0.92)',
        }
  }, [isDark])

  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas) return

    const reduced = prefersReducedMotion()
    const ctx = canvas.getContext('2d')
    if (!ctx) return

    // Pause animation when not visible to save CPU.
    const observer = new IntersectionObserver(
      (entries) => {
        runningRef.current = entries.some((e) => e.isIntersecting)
      },
      { threshold: 0.05 }
    )
    observer.observe(canvas)

    const resize = () => {
      const dpr = Math.max(1, Math.min(2, window.devicePixelRatio || 1))
      const rect = canvas.getBoundingClientRect()
      const w = Math.max(1, Math.floor(rect.width * dpr))
      const h = Math.max(1, Math.floor(rect.height * dpr))
      if (canvas.width !== w || canvas.height !== h) {
        canvas.width = w
        canvas.height = h
      }
    }

    const initParticles = () => {
      resize()
      const w = canvas.width
      const h = canvas.height

      const count = Math.max(8, Math.min(60, Math.floor(particleCount)))
      const rand = (a: number, b: number) => a + Math.random() * (b - a)

      particlesRef.current = new Array(count).fill(0).map((_, i) => {
        const baseR = rand(2.4, 4.6)
        const speed = rand(0.18, 0.52) * (Math.min(w, h) / 260)
        const angle = rand(0, Math.PI * 2)
        const type: 0 | 1 = i % 5 === 0 ? 1 : 0
        return {
          x: rand(baseR, w - baseR),
          y: rand(baseR, h - baseR),
          vx: Math.cos(angle) * speed,
          vy: Math.sin(angle) * speed,
          r: baseR,
          type,
        }
      })
    }

    initParticles()

    const step = () => {
      rafRef.current = window.requestAnimationFrame(step)
      if (reduced || !runningRef.current) return

      const w = canvas.width
      const h = canvas.height

      // No trailing: redraw clean every frame.
      ctx.clearRect(0, 0, w, h)

      // Bounding box (like the sim wireframe)
      const pad = Math.max(6, Math.min(12, Math.floor(Math.min(w, h) * 0.06)))
      ctx.strokeStyle = palette.box
      ctx.lineWidth = 1
      ctx.strokeRect(pad, pad, w - pad * 2, h - pad * 2)

      const ps = particlesRef.current

      // Lightweight "MD-ish" motion: periodic wrap + gentle repulsion.
      const boxL = w - pad * 2
      const boxH = h - pad * 2

      // Pairwise soft repulsion (small N => OK)
      for (let i = 0; i < ps.length; i++) {
        for (let j = i + 1; j < ps.length; j++) {
          const a = ps[i]
          const b = ps[j]
          let dx = b.x - a.x
          let dy = b.y - a.y

          // Minimum image in the thumbnail box
          if (dx > boxL / 2) dx -= boxL
          if (dx < -boxL / 2) dx += boxL
          if (dy > boxH / 2) dy -= boxH
          if (dy < -boxH / 2) dy += boxH

          const r0 = a.r + b.r + 6
          const r2 = dx * dx + dy * dy
          if (r2 < r0 * r0 && r2 > 1e-6) {
            const inv = 1 / Math.sqrt(r2)
            const push = (r0 - Math.sqrt(r2)) * 0.0024
            const fx = dx * inv * push
            const fy = dy * inv * push
            a.vx -= fx
            a.vy -= fy
            b.vx += fx
            b.vy += fy
          }
        }
      }

      // Keep motion lively: gentle stochastic "thermostat" kicks + speed clamp.
      const kick = 0.008 * (Math.min(w, h) / 320)
      const vMin = 0.06 * (Math.min(w, h) / 260)
      const vMax = 0.9 * (Math.min(w, h) / 260)

      for (const p of ps) {
        p.vx += (Math.random() - 0.5) * kick
        p.vy += (Math.random() - 0.5) * kick

        const v = Math.hypot(p.vx, p.vy)
        if (v > 1e-8) {
          const clamped = Math.min(vMax, Math.max(vMin, v))
          const s = clamped / v
          p.vx *= s
          p.vy *= s
        }

        p.x += p.vx
        p.y += p.vy

        // Periodic wrapping within padded box
        if (p.x < pad) p.x += boxL
        if (p.x >= pad + boxL) p.x -= boxL
        if (p.y < pad) p.y += boxH
        if (p.y >= pad + boxH) p.y -= boxH
      }

      // Draw particles (blue/red only)
      for (let i = 0; i < ps.length; i++) {
        const p = ps[i]

        ctx.beginPath()
        ctx.arc(p.x, p.y, p.r, 0, Math.PI * 2)
        ctx.fillStyle = p.type === 0 ? palette.a : palette.b
        ctx.fill()
      }
    }

    // Initial clear
    ctx.clearRect(0, 0, canvas.width, canvas.height)

    // Kick off loop
    runningRef.current = true
    rafRef.current = window.requestAnimationFrame(step)

    const onResize = () => initParticles()
    window.addEventListener('resize', onResize)

    return () => {
      window.removeEventListener('resize', onResize)
      observer.disconnect()
      if (rafRef.current) window.cancelAnimationFrame(rafRef.current)
      rafRef.current = null
    }
  }, [palette, particleCount, isDark])

  return (
    <canvas
      ref={canvasRef}
      className={className}
      style={{ backgroundColor: 'hsl(var(--card))' }}
      aria-hidden="true"
    />
  )
}
