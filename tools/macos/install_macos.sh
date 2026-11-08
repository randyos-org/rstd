cd /Volumes/macOS && mkdir -p private/tmp && cd private/tmp
curl -L -C - -o InstallAssistant.pkg 'https://swcdn.apple.com/content/downloads/13/30/140-88558-A_TN3FIE7J63/p61cqj3ktmcqilz3q64ls4ekp1os9c8qbp/InstallAssistant.pkg'
ls -d /Install*.app
cp -R "/Install macOS Sonoma.app" .
mkdir "Install macOS Sonoma.app/Contents/SharedSupport"
mv InstallAssistant.pkg "Install macOS Sonoma.app/Contents/SharedSupport/SharedSupport.dmg"
"Install macOS Sonoma.app/Contents/MacOS/InstallAssistant_springboard"
