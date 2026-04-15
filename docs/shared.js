// shared.js
document.querySelectorAll('.copy-btn').forEach(btn => {
  btn.addEventListener('click', () => {
    const pre = btn.nextElementSibling;
    const text = pre ? pre.textContent : '';
    navigator.clipboard.writeText(text).then(() => {
      btn.textContent = 'copied!';
      setTimeout(() => btn.textContent = 'copy', 1600);
    });
  });
});
