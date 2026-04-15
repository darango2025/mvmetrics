# ============================================================
# mvmetrics — Paso a paso para construir e instalar la librería
# ============================================================
# Ejecuta estos comandos en orden desde la consola R o RStudio.
# Requisito previo: tener el paquete devtools instalado.
# ============================================================

# ── PASO 0: Instalar herramientas de desarrollo ───────────────────────────────
install.packages(c("devtools", "roxygen2", "testthat", "usethis"))

# ── PASO 1: Definir la ruta de trabajo ────────────────────────────────────────
# Ajusta esta ruta a donde guardaste la carpeta mvmetrics/
pkg_path <- "~/mvmetrics"       # <-- CAMBIA ESTA RUTA

# ── PASO 2: Generar documentación con roxygen2 ────────────────────────────────
# Convierte los bloques #' en archivos .Rd dentro de man/
devtools::document(pkg_path)
# Resultado esperado: "Writing NAMESPACE" + un .Rd por función exportada

# ── PASO 3: Verificar integridad del paquete ──────────────────────────────────
devtools::check(pkg_path)
# IDEAL: 0 errors, 0 warnings, 0 notes
# (Es normal un NOTE sobre "non-standard files" la primera vez)

# ── PASO 4: Instalar el paquete localmente ────────────────────────────────────
devtools::install(pkg_path)
# Ahora puedes hacer: library(mvmetrics)

# ── PASO 5: Ejecutar los tests unitarios ──────────────────────────────────────
devtools::test(pkg_path)
# Resultado esperado: 20+ tests, todos OK / PASS

# ============================================================
# USO BÁSICO — Ejemplo completo con datos simulados
# ============================================================
library(mvmetrics)

set.seed(2026)
N <- 500; K <- 5; K_c <- 4; B <- 50

# Parámetros empíricos (Valle del Cauca calibrated)
mu_real  <- c(Tmin = 18.9, Tmax = 30.6, HR = 80.7, Rad = 3.8)
sd_real  <- c(Tmin = 2.5,  Tmax = 4.0,  HR = 8.0,  Rad = 0.8)

# Datos observados: 4 continuas + 1 binaria
obs <- cbind(
  matrix(rnorm(N * K_c, rep(mu_real, each = N), rep(sd_real, each = N)), N, K_c),
  rbinom(N, 1, 0.35)
)
colnames(obs) <- c("Tmin","Tmax","HR","Rad","Pbin")

# M1: modelo independiente (baseline)
m1_pred <- cbind(
  matrix(rnorm(N * K_c, rep(mu_real, each = N), rep(sd_real * 1.1, each = N)), N, K_c),
  runif(N, 0.3, 0.45)
)
m1_draws <- array(rnorm(N * K * B), dim = c(N, K, B))

# M3: modelo conjunto (mejor estructura de dependencia)
m3_pred <- cbind(
  matrix(rnorm(N * K_c, rep(mu_real, each = N), rep(sd_real * 0.85, each = N)), N, K_c),
  runif(N, 0.32, 0.50)
)
m3_draws <- array(rnorm(N * K * B), dim = c(N, K, B))

colnames(m1_pred) <- colnames(m3_pred) <- colnames(obs)

# ── Evaluación completa con todas las métricas ─────────────────────────────
report <- evaluate_models(
  obs        = obs,
  pred_list  = list(M1 = m1_pred, M3 = m3_pred),
  draws_list = list(M1 = m1_draws, M3 = m3_draws),
  mu_train   = mu_real,
  sd_train   = sd_real,
  K_c        = K_c,
  alpha_D    = 0.5,       # peso fijo para Métrica D
  alpha_E    = NULL,      # alpha* adaptativo para Métrica E
  run_dcor   = FALSE,     # TRUE para análisis final (más lento)
  bootstrap  = TRUE,      # CIs bootstrap
  B_boot     = 200,       # resamples (500 para análisis final)
  verbose    = TRUE
)

# ── Visualización de resultados ───────────────────────────────────────────
print(report)           # tabla completa + recomendación
summary(report)         # resumen compacto

# Gráficos (requiere ggplot2)
plot_metrics(report, type = "bars")     # barras por métrica
plot_metrics(report, type = "heatmap")  # mapa de calor de rankings
plot_metrics(report, type = "decomp")   # S_marg vs S_dep
plot_metrics(report, type = "forest")   # CIs bootstrap (forest plot)

# ── Uso de funciones individuales ─────────────────────────────────────────
# Sólo RMSE marginal
marginal_rmse(obs[, 1:K_c], m1_pred[, 1:K_c])

# Normalizar antes de B y C
obs_n    <- normalise_matrix(obs[, 1:K_c], mu_real, sd_real)

# Energy Score para M1
metric_B_energy_score(cbind(obs_n, obs[, K]), m1_draws)

# MDD para M3 con alpha sensibility
for (a in c(0.3, 0.5, 0.7)) {
  d <- metric_D_mdd(obs, m3_pred, m3_draws, K_c = K_c, alpha = a)
  cat(sprintf("alpha=%.1f  D=%.4f  S_dep=%.4f\n", a, d$total, d$S_dep))
}

# E-CVWMD: ver pesos CV y alpha*
res_E <- metric_E_cvwmd(obs, m3_pred, m3_draws,
                         mu_train = mu_real, sd_train = sd_real,
                         K_c = K_c, run_dcor = FALSE)
cat("CV weights:", round(res_E$w_cont, 3), "| binary:", round(res_E$w_bin, 3))
cat("\nalpha* =", res_E$alpha_star)

# ============================================================
# ESTRUCTURA DE ARCHIVOS DEL PAQUETE
# ============================================================
# mvmetrics/
# ├── DESCRIPTION          — metadatos del paquete
# ├── NAMESPACE            — funciones exportadas (generado por roxygen2)
# ├── LICENSE
# ├── R/
# │   ├── utils.R          — normalise_matrix + helpers internos
# │   ├── marginal_metrics.R — marginal_rmse, marginal_mae, marginal_logloss,
# │   │                         marginal_accuracy, metric_A
# │   ├── joint_metrics.R  — metric_B_energy_score, metric_C_variogram_score,
# │   │                         metric_D_mdd, metric_E_cvwmd
# │   ├── evaluate_models.R — evaluate_models (wrapper principal),
# │   │                         print/summary methods
# │   ├── bootstrap.R      — bootstrap_metrics
# │   ├── recommend.R      — recommend_metric
# │   └── plot_metrics.R   — plot_metrics (requires ggplot2)
# └── tests/
#     └── testthat/
#         └── test-metrics.R — 20+ unit tests
