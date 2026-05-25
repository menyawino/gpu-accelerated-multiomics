#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$ROOT_DIR/results/anydesk/connect"
mkdir -p "$OUT_DIR"

wsl() {
  /mnt/c/Windows/System32/wsl.exe -d "$1" --cd /root -- bash -lc "$2"
}

new_session() {
  local distro="$1"
  local sock="$2"

  wsl "$distro" "pkill -f '^tmate' || true"
  wsl "$distro" "rm -f '$sock'"
  wsl "$distro" "tmate -S '$sock' new-session -d -c \"\$HOME\""
  wsl "$distro" "tmate -S '$sock' wait tmate-ready"

  local rw ro web
  rw="$(wsl "$distro" "tmate -S '$sock' display -p '#{tmate_ssh}'")"
  ro="$(wsl "$distro" "tmate -S '$sock' display -p '#{tmate_ssh_ro}'")"
  web="$(wsl "$distro" "tmate -S '$sock' display -p '#{tmate_web}'")"

  printf '%s\n%s\n%s\n' "$rw" "$ro" "$web"
}

readarray -t EYAD_LINES < <(new_session "Ubuntu" "/tmp/tmate-eyad.sock")
readarray -t MAHMOUD_LINES < <(new_session "Ubuntu-clone" "/tmp/tmate-mahmoud.sock")

EYAD_RW="${EYAD_LINES[0]}"
EYAD_RO="${EYAD_LINES[1]}"
EYAD_WEB="${EYAD_LINES[2]}"
MAHMOUD_RW="${MAHMOUD_LINES[0]}"
MAHMOUD_RO="${MAHMOUD_LINES[1]}"
MAHMOUD_WEB="${MAHMOUD_LINES[2]}"

EYAD_HOST="${EYAD_RW#ssh }"
EYAD_HOST="${EYAD_HOST#\"}"
EYAD_HOST="${EYAD_HOST%\"}"
MAHMOUD_HOST="${MAHMOUD_RW#ssh }"
MAHMOUD_HOST="${MAHMOUD_HOST#\"}"
MAHMOUD_HOST="${MAHMOUD_HOST%\"}"

cat > "$OUT_DIR/eyad_connect.sh" <<EOF
#!/usr/bin/env bash
exec ssh $EYAD_HOST
EOF

cat > "$OUT_DIR/mahmoud_connect.sh" <<EOF
#!/usr/bin/env bash
exec ssh $MAHMOUD_HOST
EOF

chmod +x "$OUT_DIR/eyad_connect.sh" "$OUT_DIR/mahmoud_connect.sh"

cat > "$OUT_DIR/current_links.txt" <<EOF
EYAD_RW=$EYAD_RW
EYAD_RO=$EYAD_RO
EYAD_WEB=$EYAD_WEB
MAHMOUD_RW=$MAHMOUD_RW
MAHMOUD_RO=$MAHMOUD_RO
MAHMOUD_WEB=$MAHMOUD_WEB
EOF

ORIGIN_URL="$(git -C "$ROOT_DIR" remote get-url origin 2>/dev/null || true)"
RAW_BASE=""
if [[ "$ORIGIN_URL" =~ ^https://github.com/([^/]+)/([^/.]+)(\.git)?$ ]]; then
  RAW_BASE="https://raw.githubusercontent.com/${BASH_REMATCH[1]}/${BASH_REMATCH[2]}/main"
fi

if [[ -n "$RAW_BASE" ]]; then
  cat > "$OUT_DIR/one_command_connect.txt" <<EOF
Eyad one-command connect:
bash <(curl -fsSL $RAW_BASE/results/anydesk/connect/eyad_connect.sh)

Mahmoud one-command connect:
bash <(curl -fsSL $RAW_BASE/results/anydesk/connect/mahmoud_connect.sh)
EOF
fi

echo "Rotation complete."
echo "$EYAD_RW"
echo "$MAHMOUD_RW"
echo "Artifacts: $OUT_DIR"
