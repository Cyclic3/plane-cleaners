#!/bin/sh
set -e

if [ "$#" -ne 2 ]
then
  echo "Usage: $0 <binary> <destination>"
  exit 1
fi

export DEPLOY_DIR="$(dirname $(readlink -f "${0}"))"
export SOURCE_BIN="$(readlink -f "$1")"
export OUTPUT="$(readlink -f "$2")"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf -- "$WORK_DIR"' EXIT

cd "$WORK_DIR"
mkdir AppDir
cd AppDir

export APPRUN_URL="https://github.com/AppImage/AppImageKit/releases/download/continuous/AppRun-x86_64"
export LD_SO_PATH="/lib64/ld-linux-x86-64.so.2"

curl -L "$APPRUN_URL" -o AppRun
chmod +x AppRun

cat << EOF > plane-cleaners.desktop
[Desktop Entry]
Encoding=UTF-8
Version=1.0
Type=Application
Terminal=true
Exec=start.sh
Name=plane-cleaners
Icon=empty
Categories=Utility;
EOF
cp "$DEPLOY_DIR/empty.png" .

mkdir -p usr/bin
cd usr/bin

cat << EOF > start.sh
#!/bin/sh
export LIB_PATH="\$(dirname \$(readlink -f "\${0}"))"
"\$LIB_PATH/ld.so" --library-path "\$LIB_PATH" "\$LIB_PATH/$(basename $1)" "\$@"
EOF

chmod +x start.sh

cp "$SOURCE_BIN" .

export FILES="$(ldd $1 | awk '{if (NF >= 3) print $3}')"

ldd "$SOURCE_BIN" | awk '{if (NF >= 3) print $3}' | while read file
do
        cp $file .
done

cp "$LD_SO_PATH" ld.so

cd "$WORK_DIR"
wget https://github.com/AppImage/AppImageKit/releases/download/continuous/appimagetool-x86_64.AppImage
chmod +x appimagetool-x86_64.AppImage
./appimagetool-x86_64.AppImage AppDir "$OUTPUT"
