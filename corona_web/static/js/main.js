// Wait for DOM to fully load before running scripts
document.addEventListener("DOMContentLoaded", function () {
    const form = document.getElementById("searchForm");
    const typeSelect = document.getElementById("typeSelect");

    // Check if the form exists on the page
    if (form) {
        // Attach submit event listener to the form
        form.addEventListener("submit", function (e) {
            e.preventDefault();  // Prevent default form submission
            const type = typeSelect.value;  // Get selected data type (genes/proteins/pdbs)
            const virus = document.getElementById("virusSelect").value;  // Get selected virus
            // Redirect to the appropriate page with virus query parameter
            window.location.href = `/${type}?virus=${encodeURIComponent(virus)}`;
        });
    }
});