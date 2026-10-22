// Render an HTML file to PNG with the installed Microsoft Edge (needs: npm i playwright-core).
// usage: node tools/render_png.js <input.html> <output.png> [width=1920] [height=1080] [scale=2]
const { chromium } = require("playwright-core");
const path = require("path");

(async () => {
  const [inp, out, w = "1920", h = "1080", s = "2"] = process.argv.slice(2);
  const browser = await chromium.launch({ channel: "msedge", headless: true });
  const page = await (await browser.newContext({ viewport: { width: +w, height: +h }, deviceScaleFactor: +s })).newPage();
  await page.goto("file:///" + path.resolve(inp).replace(/\\/g, "/"), { waitUntil: "load" });
  await page.evaluate(() => document.fonts.ready);
  await page.waitForTimeout(500);
  const overflow = await page.evaluate(() => ({ w: document.documentElement.scrollWidth, h: document.documentElement.scrollHeight }));
  await page.screenshot({ path: out });
  console.log(`saved ${out} (content ${overflow.w}x${overflow.h} css px)`);
  await browser.close();
})().catch((e) => { console.error(e); process.exit(1); });
