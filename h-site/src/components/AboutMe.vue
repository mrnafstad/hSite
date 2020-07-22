<template>
	<div id="about">
		<li v-for="(element, idx) in Om" v-bind:key="idx"> {{ element.desc }} </li>
	</div>
</template>

<script>
	module.exports = {
		data: function() {
			return {
				Om: [
				{
					id: 1,
					desc: "I'm a 29 year old developer. I got my MSc in astrophysics in the spring of 2019. During post graduation job hunting I branched out and discoverred a passion for web development, particularly frontend."
				},
				{ 
					id: 2,
					desc: "Therefore I've spent the last year back in the classroom taking courses to become a better over-all programmer."
				},
				{
					id: 3,
					desc: "Some of the courses I've taken covered topics like functional programming, Java, SQL, NoSQL (MongoDB), OpenMP and MPI. On the side I've tried to learn the Vue framework, which I'm building this website with."
				}],
			};
		}
	};
</script>

<style>
	li {
		padding: 5px;
		list-style-type: none;
	}
	#about {
		width: 2;
		margin-left: auto;
		margin-right: auto;
	}
</style>