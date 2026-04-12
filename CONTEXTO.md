---
name: ahorro-tokens
description: Reglas para reducir tokens en Claude Code
---

#
meta:
  title: "Reglas para Claude Code — Ahorra Tokens"
  instruction: "Copia este contenido en `CLAUDE.md` en la raiz de tu proyecto o en `~/.claude/CLAUDE.md` para que aplique a todos tus proyectos."
rules:
  - number: 1
    title: "No programar sin contexto"
    description: |
      ANTES de escribir codigo: lee los archivos relevantes, revisa git log, entiende la arquitectura.
      Si no tienes contexto suficiente, pregunta. No asumas.
  - number: 2
    title: "Respuestas cortas"
    description: |
      Responde en 1-3 oraciones. Sin preambulos, sin resumen final.
      No repitas lo que el usuario dijo. No expliques lo obvio.
      Codigo habla por si mismo: no narres cada linea que escribes.
  - number: 3
    title: "No reescribir archivos completos"
    description: |
      Usa Edit (reemplazo parcial), NUNCA Write para archivos existentes salvo que el cambio sea >80% del archivo.
      Cambia solo lo necesario. No "limpies" codigo alrededor del cambio.
  - number: 4
    title: "No releer archivos ya leidos"
    description: |
      Si ya leiste un archivo en esta conversacion, no lo vuelvas a leer salvo que haya cambiado.
      Toma notas mentales de lo importante en tu primera lectura.
  - number: 5
    title: "Validar antes de declarar hecho"
    description: |
      Despues de un cambio: compila, corre tests, o verifica que funciona.
      Nunca digas "listo" sin evidencia de que funciona.
  - number: 6
    title: "Cero charla aduladora"
    description: |
      No digas "Excelente pregunta", "Gran idea", "Perfecto", etc.
      No halagues al usuario. Ve directo al trabajo.
  - number: 7
    title: "Soluciones simples"
    description: |
      Implementa lo minimo que resuelve el problema. Nada mas.
      No agregues abstracciones, helpers, tipos, validaciones, ni features que no se pidieron.
      3 lineas repetidas > 1 abstraccion prematura.
  - number: 8
    title: "No pelear con el usuario"
    description: |
      Si el usuario dice "hazlo asi", hazlo asi. No debatas salvo riesgo real de seguridad o perdida de datos.
      Si discrepas, menciona tu concern en 1 oracion y procede con lo que pidio.
  - number: 9
    title: "Leer solo lo necesario"
    description: |
      No leas archivos completos si solo necesitas una seccion. Usa offset y limit.
      Si sabes la ruta exacta, usa Read directo. No hagas Glob + Grep + Read cuando Read basta.
  - number: 10
    title: "No narrar el plan antes de ejecutar"
    description: |
      No digas "Voy a leer el archivo, luego modificar la funcion, luego compilar...". Solo hazlo.
      El usuario ve tus tool calls. No necesita un preview en texto.
  - number: 11
    title: "Paralelizar tool calls"
    description: |
      Si necesitas leer 3 archivos independientes, lee los 3 en un solo mensaje, no uno por uno.
      Menos roundtrips = menos tokens de contexto acumulado.
  - number: 12
    title: "No duplicar codigo en la respuesta"
    description: |
      Si ya editaste un archivo, no copies el resultado en tu respuesta. El usuario lo ve en el diff.
      Si creaste un archivo, no lo muestres entero en texto tambien.
  - number: 13
    title: "No usar Agent cuando Grep/Read basta"
    description: |
      Agent duplica todo el contexto en un subproceso. Solo usalo para busquedas amplias o tareas complejas.
      Para buscar una funcion o archivo especifico, usa Grep o Glob directo.
