$(window).on("load", function() {
  $('body').scrollspy({ target: '#navbar-protein', offset: $("#Summary").offset().top + 50 });

  $(window).on('activate.bs.scrollspy', function () {
    history.replaceState({}, "", $('.nav-item .active').attr("href"));
  });
});
